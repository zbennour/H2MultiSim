"""
Binary H2 + cushion-gas composition support shared by the 2D and 3D tiers.

The grid models carry one H2 mole fraction per cell.  Cell properties come
from the single-gas EOS when the run is pure H2 (unchanged fast path) and from
a :class:`PengRobinsonMixture` otherwise; composition is transported with an
implicit upwind advection–dispersion scheme in advective form (as in
``UHSRadial1D``), which is exactly H2-conservative at a converged Picard step.
"""

from __future__ import annotations

from h2multisim.thermo import (
    GAS_COMPONENTS,
    M_H2,
    PengRobinsonMixture,
    isentropic_exponent,
    viscosity_h2,
)


def validate_cushion_props(props: object, where: str) -> None:
    """Normalise ``props.cushion_gas`` and check ``initial_h2_mole_fraction``."""
    cg = str(props.cushion_gas).strip().upper()
    if cg not in GAS_COMPONENTS:
        raise ValueError(
            f"{where}: unknown cushion_gas {props.cushion_gas!r}; choose from {list(GAS_COMPONENTS)}"
        )
    props.cushion_gas = cg
    y0 = props.initial_h2_mole_fraction
    if y0 is not None and not (0.0 <= y0 <= 1.0):
        raise ValueError(f"{where}: initial_h2_mole_fraction must be in [0, 1]")


class BinaryCompositionMixin:
    """Composition-aware cell properties; requires ``self.eos`` and ``self.props``."""

    def _init_composition(self) -> None:
        props = self.props
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
        self._m_cg: float = GAS_COMPONENTS[cg].molar_mass
        self._transport_failures: int = 0

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

    def _w_from_y(self, y: float) -> float:
        """H2 mass fraction from H2 mole fraction."""
        return y * M_H2 / max(y * M_H2 + (1.0 - y) * self._m_cg, 1e-30)

    def _y_from_w(self, w: float) -> float:
        """H2 mole fraction from H2 mass fraction."""
        n_h2 = w / M_H2
        n_cg = (1.0 - w) / self._m_cg
        return n_h2 / max(n_h2 + n_cg, 1e-30)


def add_transport_link(
    coup: list[list[tuple[int, float]]], c1: int, c2: int, F: float, G: float
) -> None:
    """
    Register an interface between flat cells ``c1`` and ``c2``.  ``F`` [kg/s] is
    the mass flow from c1 to c2 and ``G`` [kg/s] the dispersive conductance; each
    cell couples to its neighbour through upwind inflow plus dispersion.
    """
    coup[c2].append((c1, max(F, 0.0) + G))
    coup[c1].append((c2, max(-F, 0.0) + G))


def solve_binary_transport(
    m_old: list[float],
    w_old: list[float],
    coup: list[list[tuple[int, float]]],
    q_inj: list[float],
    dt: float,
    max_sweeps: int = 500,
    tol: float = 1e-12,
) -> tuple[list[float], bool]:
    """
    Implicit advective-form transport of the H2 mass fraction w on flat cells:

        (m_old/dt + Σ c_nb + q_inj)·w − Σ c_nb·w_nb = m_old/dt·w_old + q_inj·1

    (injected gas is pure H2; production removes the cell composition and needs
    no term).  The matrix is a diagonally dominant M-matrix, so Gauss–Seidel
    converges.  Returns ``(w_new, converged)``.
    """
    n = len(w_old)
    w = list(w_old)
    base = [m_old[c] / dt for c in range(n)]
    diag = [base[c] + q_inj[c] + sum(cf for _, cf in coup[c]) for c in range(n)]
    for _ in range(max_sweeps):
        max_dw = 0.0
        for c in range(n):
            rhs = base[c] * w_old[c] + q_inj[c]
            for nb, cf in coup[c]:
                rhs += cf * w[nb]
            w_val = min(max(rhs / diag[c], 0.0), 1.0)
            d = abs(w_val - w[c])
            if d > max_dw:
                max_dw = d
            w[c] = w_val
        if max_dw < tol:
            return w, True
    return w, False
