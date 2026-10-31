"""Configuration loading for H2MultiSim.

Allows a run to be fully described by a JSON (or YAML, if PyYAML is
installed) file instead of hard-coded scenario functions.  Every field maps
directly onto the dataclasses in :mod:`h2multisim.models`, so the file is
the single source of truth for reservoir properties, wells and schedule.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

from h2multisim.models.one_d import (
    RadialReservoirProps,
    WellControl,
)
from h2multisim.models.three_d import (
    Reservoir3DProps,
    Schedule3D,
    WellSpec3D,
)
from h2multisim.models.two_d import (
    Reservoir2DProps,
    Schedule2D,
    WellSpec2D,
)
from h2multisim.models.zero_d import (
    ReservoirProperties,
    SchedulePeriod,
)
from h2multisim.thermo import (
    IsothermModel,
    PoreVolumeModel,
    RelativePermeabilityModel,
)
from h2multisim.units import mD_to_m2


def _load_raw(path: str | os.PathLike) -> dict[str, Any]:
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"config file not found: {p}")
    text = p.read_text(encoding="utf-8")
    if p.suffix.lower() in (".yaml", ".yml"):
        try:
            import yaml  # type: ignore
        except ImportError as exc:  # pragma: no cover
            raise RuntimeError("YAML support requires PyYAML: pip install pyyaml") from exc
        data = yaml.safe_load(text)
    else:
        data = json.loads(text)
    if not isinstance(data, dict):
        raise ValueError("config must be a JSON/YAML object")
    return data


_PRESSURE_UNITS: tuple[tuple[str, float], ...] = (
    ("bara", 1.0e5),
    ("bar", 1.0e5),
    ("mpa", 1.0e6),
    ("kpa", 1.0e3),
    ("psia", 6894.757293168),
    ("psi", 6894.757293168),
    ("pa", 1.0),
)


def _pa(v: Any) -> float:
    if isinstance(v, str):
        s = v.strip().lower()
        for suffix, factor in _PRESSURE_UNITS:
            if s.endswith(suffix):
                return float(s[: -len(suffix)].strip()) * factor
        raise ValueError(f"cannot parse pressure: {v!r} (units: bar, MPa, kPa, Pa, psi)")
    return float(v)


def _perm(d: dict[str, Any]) -> float:
    if "k" in d:
        return float(d["k"])
    if "permeability" in d:
        return float(d["permeability"])
    if "permeability_mD" in d:
        return mD_to_m2(float(d["permeability_mD"]))
    return 0.0


def _rel_perm(d: dict[str, Any]) -> RelativePermeabilityModel | None:
    rp_raw = d.get("rel_perm")
    if isinstance(rp_raw, dict):
        return RelativePermeabilityModel(
            model=str(rp_raw.get("model", "corey")),
            swr=float(rp_raw.get("swr", d.get("swr", 0.0))),
            sgc=float(rp_raw.get("sgc", d.get("sgc", 0.0))),
            krg_max=float(rp_raw["krg_max"]) if rp_raw.get("krg_max") is not None else None,
            krw_max=float(rp_raw.get("krw_max", 1.0)),
            n_g=float(rp_raw.get("n_g", 2.0)),
            n_w=float(rp_raw.get("n_w", 3.0)),
            brooks_corey_lambda=float(rp_raw.get("brooks_corey_lambda", 2.0)),
            entry_pressure_Pa=_pa(rp_raw.get("entry_pressure_Pa", 5.0e3)),
            active_gas_saturation=(
                float(rp_raw["active_gas_saturation"])
                if rp_raw.get("active_gas_saturation") is not None
                else None
            ),
        )
    return None


def _props_0d(d: dict[str, Any]) -> ReservoirProperties:
    area = float(d["area"])
    thickness = float(d["thickness"])
    porosity = float(d["porosity"])
    ct = float(d["total_compressibility"])
    p_ref = _pa(d["reference_pressure"])
    rel_perm = _rel_perm(d)
    # Default pore-model volume is the gas-accessible pore volume, matching
    # UHSReservoir0D when no pore model is given.
    sat_factor = (
        rel_perm
        or RelativePermeabilityModel(swr=float(d.get("swr", 0.0)), sgc=float(d.get("sgc", 0.0)))
    ).effective_porosity_factor()
    v0_default = area * thickness * porosity * sat_factor

    pm_raw = d.get("pore_model")
    pore_model: PoreVolumeModel | None = None
    if isinstance(pm_raw, dict):
        pore_model = PoreVolumeModel(
            model=str(pm_raw.get("model", "elastic")),
            V0=float(pm_raw.get("V0", v0_default)),
            P0=_pa(pm_raw.get("P0", p_ref)),
            total_compressibility=float(pm_raw.get("total_compressibility", ct)),
            P_gate=_pa(pm_raw.get("P_gate", 5e6)),
            k_gate=float(pm_raw.get("k_gate", 1e-6)),
            V_ratio=float(pm_raw.get("V_ratio", 2.0)),
        )
    elif isinstance(pm_raw, str):
        pore_model = PoreVolumeModel(
            model=pm_raw,
            V0=v0_default,
            P0=p_ref,
            total_compressibility=ct,
        )

    iso_raw = d.get("isotherm")
    isotherm: IsothermModel | None = None
    if isinstance(iso_raw, dict):
        isotherm = IsothermModel(
            model=str(iso_raw.get("model", "langmuir")),
            n_max1=float(iso_raw.get("n_max1", 20.0)),
            K0_1=float(iso_raw.get("K0_1", 1e-7)),
            n_max2=float(iso_raw.get("n_max2", 5.0)),
            K0_2=float(iso_raw.get("K0_2", 5e-8)),
            T_ref=float(iso_raw.get("T_ref", 298.0)),
            dH1=float(iso_raw.get("dH1", -6000.0)),
            dH2=float(iso_raw.get("dH2", -4000.0)),
        )
    elif isinstance(iso_raw, str):
        isotherm = IsothermModel(model=iso_raw)

    return ReservoirProperties(
        area=area,
        thickness=thickness,
        porosity=porosity,
        total_compressibility=ct,
        permeability=_perm(d),
        reference_pressure=p_ref,
        reference_temperature=float(d["reference_temperature"]),
        well_radius=float(d["well_radius"]),
        drainage_radius=float(d["drainage_radius"]),
        is_cavern=bool(d.get("is_cavern", False)),
        n_wells=int(d.get("n_wells", 1)),
        pore_model=pore_model,
        isotherm=isotherm,
        m_sorbent=float(d.get("m_sorbent", 0.0)),
        volume_translation=float(
            d.get("volume_translation", -4.0e-6 if d.get("peneloux", False) else 0.0)
        ),
        peneloux=bool(d.get("peneloux", False)),
        cushion_gas=str(d.get("cushion_gas", "H2")),
        initial_composition=d.get("initial_composition"),
        injected_composition=d.get("injected_composition"),
        plume_segregation=float(d.get("plume_segregation", 0.0)),
        swr=float(d.get("swr", 0.0)),
        sgc=float(d.get("sgc", 0.0)),
        rel_perm=rel_perm,
        thermal_model=str(d.get("thermal_model", "energy_balance")),
        supply_pressure=_pa(d.get("supply_pressure", 30.0e5)),
        compressor_efficiency=float(d.get("compressor_efficiency", 0.75)),
        compressor_max_stage_ratio=float(d.get("compressor_max_stage_ratio", 3.0)),
        rock_heat_capacity=(
            float(d["rock_heat_capacity"]) if d.get("rock_heat_capacity") is not None else None
        ),
        injection_temperature=(
            float(d["injection_temperature"])
            if d.get("injection_temperature") is not None
            else None
        ),
    )


def _props_1d(d: dict[str, Any]) -> RadialReservoirProps:
    return RadialReservoirProps(
        r_well=float(d["r_well"]),
        r_ext=float(d["r_ext"]),
        h=float(d["h"]),
        phi=float(d["phi"]),
        k=_perm(d),
        Ct=float(d["Ct"]),
        P_init=_pa(d["P_init"]),
        T_init=float(d["T_init"]),
        skin=float(d.get("skin", 0.0)),
        outer_bc=str(d.get("outer_bc", "no_flow")),
        tau_days=float(d.get("tau_days", 30.0)),
        volume_translation=float(
            d.get("volume_translation", -4.0e-6 if d.get("peneloux", False) else 0.0)
        ),
        peneloux=bool(d.get("peneloux", False)),
        cushion_gas=str(d.get("cushion_gas", "H2")),
        initial_h2_mole_fraction=(
            float(d["initial_h2_mole_fraction"])
            if d.get("initial_h2_mole_fraction") is not None
            else None
        ),
        dispersivity_m=float(d.get("dispersivity_m", 0.5)),
        molecular_diffusivity_m2_s=float(d.get("molecular_diffusivity_m2_s", 1.0e-5)),
        swr=float(d.get("swr", 0.0)),
        sgc=float(d.get("sgc", 0.0)),
        rel_perm=_rel_perm(d),
        rock_heat_capacity=float(d.get("rock_heat_capacity", 2.2e6)),
    )


def _props_2d(d: dict[str, Any]) -> Reservoir2DProps:
    k_field: list[list[float]] | None = None
    if d.get("k_field") is not None:
        k_field = [[float(val) for val in row] for row in d["k_field"]]
    elif d.get("k_field_mD") is not None:
        k_field = [[mD_to_m2(float(val)) for val in row] for row in d["k_field_mD"]]

    return Reservoir2DProps(
        Lx=float(d["Lx"]),
        Ly=float(d["Ly"]),
        h=float(d["h"]),
        phi=float(d["phi"]),
        k=_perm(d),
        Ct=float(d["Ct"]),
        P_init=_pa(d["P_init"]),
        T_init=float(d["T_init"]),
        nx=int(d["nx"]),
        ny=int(d["ny"]),
        tau_days=float(d.get("tau_days", 30.0)),
        k_field=k_field,
        volume_translation=float(
            d.get("volume_translation", -4.0e-6 if d.get("peneloux", False) else 0.0)
        ),
        peneloux=bool(d.get("peneloux", False)),
        swr=float(d.get("swr", 0.0)),
        sgc=float(d.get("sgc", 0.0)),
        rel_perm=_rel_perm(d),
        rock_heat_capacity=float(d.get("rock_heat_capacity", 2.2e6)),
        cushion_gas=str(d.get("cushion_gas", "H2")),
        initial_h2_mole_fraction=(
            float(d["initial_h2_mole_fraction"])
            if d.get("initial_h2_mole_fraction") is not None
            else None
        ),
        dispersivity_m=float(d.get("dispersivity_m", 0.5)),
        molecular_diffusivity_m2_s=float(d.get("molecular_diffusivity_m2_s", 1.0e-5)),
    )


def _props_3d(d: dict[str, Any]) -> Reservoir3DProps:
    k_layers: list[float] | None = None
    if d.get("k_layers") is not None:
        k_layers = [float(val) for val in d["k_layers"]]
    elif d.get("k_layers_mD") is not None:
        k_layers = [mD_to_m2(float(val)) for val in d["k_layers_mD"]]

    return Reservoir3DProps(
        Lx=float(d["Lx"]),
        Ly=float(d["Ly"]),
        Lz=float(d["Lz"]),
        phi=float(d["phi"]),
        k=_perm(d),
        Ct=float(d["Ct"]),
        P_init=_pa(d["P_init"]),
        T_init=float(d["T_init"]),
        nx=int(d["nx"]),
        ny=int(d["ny"]),
        nz=int(d["nz"]),
        tau_days=float(d.get("tau_days", 30.0)),
        k_layers=k_layers,
        volume_translation=float(
            d.get("volume_translation", -4.0e-6 if d.get("peneloux", False) else 0.0)
        ),
        peneloux=bool(d.get("peneloux", False)),
        swr=float(d.get("swr", 0.0)),
        sgc=float(d.get("sgc", 0.0)),
        rel_perm=_rel_perm(d),
        rock_heat_capacity=float(d.get("rock_heat_capacity", 2.2e6)),
        cushion_gas=str(d.get("cushion_gas", "H2")),
        initial_h2_mole_fraction=(
            float(d["initial_h2_mole_fraction"])
            if d.get("initial_h2_mole_fraction") is not None
            else None
        ),
        dispersivity_m=float(d.get("dispersivity_m", 0.5)),
        molecular_diffusivity_m2_s=float(d.get("molecular_diffusivity_m2_s", 1.0e-5)),
    )


def _sched_0d(items: list[dict[str, Any]]) -> list[SchedulePeriod]:
    out: list[SchedulePeriod] = []
    for it in items:
        out.append(
            SchedulePeriod(
                name=str(it["name"]),
                mode=str(it["mode"]),
                duration_days=float(it["duration_days"]),
                bhp_target=_pa(it["bhp_target"]),
                max_mass_rate=float(it["max_mass_rate"]),
            )
        )
    return out


def _sched_2d(items: list[dict[str, Any]]) -> list[Schedule2D]:
    out: list[Schedule2D] = []
    for it in items:
        out.append(
            Schedule2D(
                name=str(it["name"]),
                mode=str(it["mode"]),
                duration_d=float(it["duration_days"]),
                bhp_target=_pa(it["bhp_target"]),
                q_max=float(it["q_max"]),
                well_modes=it.get("well_modes"),
            )
        )
    return out


def _sched_3d(items: list[dict[str, Any]]) -> list[Schedule3D]:
    out: list[Schedule3D] = []
    for it in items:
        out.append(
            Schedule3D(
                name=str(it["name"]),
                mode=str(it["mode"]),
                duration_d=float(it["duration_days"]),
                bhp_target=_pa(it["bhp_target"]),
                q_max=float(it["q_max"]),
                well_modes=it.get("well_modes"),
            )
        )
    return out


def load_config(path: str | os.PathLike) -> dict[str, Any]:
    cfg = _load_raw(path)
    if "scale" not in cfg:
        raise ValueError("config requires 'scale' (0d/1d/2d/3d)")
    cfg["scale"] = str(cfg["scale"]).lower()
    if cfg["scale"] not in ("0d", "1d", "2d", "3d"):
        raise ValueError(f"unknown scale: {cfg['scale']}")
    return cfg


def build_from_config(cfg: dict[str, Any]) -> tuple[str, Any, Any, Any, float, str]:
    scale = cfg["scale"]
    out_dir = str(cfg.get("out_dir", "output"))
    dt_days = float(cfg.get("dt_days", 1.0))

    if scale == "0d":
        props: Any = _props_0d(cfg["reservoir"])
        sched: Any = _sched_0d(cfg["schedule"])
        return scale, props, None, sched, dt_days, out_dir

    if scale == "1d":
        props = _props_1d(cfg["reservoir"])
        sched = []
        for it in cfg["schedule"]:
            sched.append(
                WellControl(
                    name=str(it.get("name", "W")),
                    mode=str(it["mode"]),
                    duration_d=float(it["duration_days"]),
                    bhp_target=_pa(it["bhp_target"]),
                    q_max=float(it["q_max"]),
                )
            )
        return scale, props, None, sched, dt_days, out_dir

    if scale == "2d":
        props = _props_2d(cfg["reservoir"])
        wells: list[WellSpec2D] = []
        for w in cfg.get("wells", []):
            wells.append(
                WellSpec2D(
                    i=int(w["i"]),
                    j=int(w["j"]),
                    r_w=float(w.get("r_w", 0.1)),
                    name=str(w.get("name", "W")),
                )
            )
        sched2: Any = _sched_2d(cfg["schedule"])
        return scale, props, wells, sched2, dt_days, out_dir

    if scale == "3d":
        props = _props_3d(cfg["reservoir"])
        wells3: list[WellSpec3D] = []
        for w in cfg.get("wells", []):
            wells3.append(
                WellSpec3D(
                    i=int(w["i"]),
                    j=int(w["j"]),
                    k_min=int(w["k_min"]),
                    k_max=int(w["k_max"]),
                    r_w=float(w.get("r_w", 0.1)),
                    name=str(w.get("name", "W")),
                )
            )
        sched3: Any = _sched_3d(cfg["schedule"])
        return scale, props, wells3, sched3, dt_days, out_dir

    raise ValueError(f"unknown scale: {scale}")
