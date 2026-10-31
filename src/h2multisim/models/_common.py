"""Shared helpers for the h2multisim solver tiers."""

from __future__ import annotations

WELL_MODES: tuple[str, ...] = ("inject", "withdraw", "shut_in")


def validate_mode(mode: str, where: str) -> str:
    """Return ``mode`` if it is a known well mode, otherwise raise ``ValueError``."""
    if mode not in WELL_MODES:
        raise ValueError(
            f"{where}: unknown well mode {mode!r}. Choose from {', '.join(WELL_MODES)}"
        )
    return mode
