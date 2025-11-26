"""
Topological phase-locking energetics.

Implements a Josephson-like energy:

    E(α) = E0 [1 - cos(n1 n2 α)]

and associated stiffness.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class LockingParameters:
    """Parameters for the phase-locking sector."""

    n1: int
    n2: int
    E0: float

    @property
    def nprod(self) -> int:
        """Product of the two winding numbers n1 n2."""
        return self.n1 * self.n2

    @property
    def stiffness(self) -> float:
        """Effective small-angle stiffness K_eff = E0 (n1 n2)^2."""
        return self.E0 * (self.nprod**2)


def energy(alpha: np.ndarray, p: LockingParameters) -> np.ndarray:
    """Full energy E(α) = E0 [1 - cos(n1 n2 α)]."""
    alpha = np.asarray(alpha, dtype=float)
    return p.E0 * (1.0 - np.cos(p.nprod * alpha))


def torque(alpha: np.ndarray, p: LockingParameters) -> np.ndarray:
    """Derivative dE/dα = E0 n1 n2 sin(n1 n2 α)."""
    alpha = np.asarray(alpha, dtype=float)
    return p.E0 * p.nprod * np.sin(p.nprod * alpha)


def small_angle_energy(alpha: np.ndarray, p: LockingParameters) -> np.ndarray:
    """Quadratic approximation near α = 0: (1/2) K_eff α²."""
    alpha = np.asarray(alpha, dtype=float)
    return 0.5 * p.stiffness * alpha**2
