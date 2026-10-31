"""Shared helpers for the h2multisim solver tiers."""

from __future__ import annotations

import math

WELL_MODES: tuple[str, ...] = ("inject", "withdraw", "shut_in")


def validate_mode(mode: str, where: str) -> str:
    """Return ``mode`` if it is a known well mode, otherwise raise ``ValueError``."""
    if mode not in WELL_MODES:
        raise ValueError(
            f"{where}: unknown well mode {mode!r}. Choose from {', '.join(WELL_MODES)}"
        )
    return mode


def check_inputs(
    where: str,
    *,
    positive: dict[str, float] | None = None,
    non_negative: dict[str, float] | None = None,
    fraction: dict[str, float] | None = None,
) -> None:
    """
    Validate physical inputs, raising ``ValueError`` with the offending field.

    ``positive``: must be finite and > 0; ``non_negative``: finite and >= 0;
    ``fraction``: finite and in (0, 1].
    """
    for name, val in (positive or {}).items():
        if not (math.isfinite(val) and val > 0.0):
            raise ValueError(f"{where}: {name} must be > 0, got {val!r}")
    for name, val in (non_negative or {}).items():
        if not (math.isfinite(val) and val >= 0.0):
            raise ValueError(f"{where}: {name} must be >= 0, got {val!r}")
    for name, val in (fraction or {}).items():
        if not (math.isfinite(val) and 0.0 < val <= 1.0):
            raise ValueError(f"{where}: {name} must be in (0, 1], got {val!r}")
