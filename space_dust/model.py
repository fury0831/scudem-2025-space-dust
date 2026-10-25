# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Fengrui Zhang, Yiye Huang, Zening Wang
"""Units and diagnostics around the original constant-density model."""
import numpy as np

from .historical import C_D, G, L_V, LAMBDA, RHO_P, particle_model


def mass_from_radius(radius_m):
    radius = np.asarray(radius_m, dtype=float)
    if np.any(~np.isfinite(radius)) or np.any(radius <= 0):
        raise ValueError("Radius must be positive and finite")
    return (4.0 / 3.0) * np.pi * RHO_P * radius**3


def radius_from_mass(mass_kg):
    mass = np.asarray(mass_kg, dtype=float)
    if np.any(~np.isfinite(mass)) or np.any(mass <= 0):
        raise ValueError("Mass must be positive and finite")
    return np.cbrt(3.0 * mass / (4.0 * np.pi * RHO_P))


def initial_state(radius_mm, speed_kms, altitude_km=100.0):
    values = np.asarray([radius_mm, speed_kms, altitude_km], dtype=float)
    if not np.isfinite(values).all() or radius_mm <= 0 or speed_kms < 0 or altitude_km < 0:
        raise ValueError("Use finite positive radius and nonnegative speed/altitude")
    return np.array([altitude_km * 1e3, speed_kms * 1e3,
                     float(mass_from_radius(radius_mm * 1e-3))])


def positive_rhs(t, state):
    """Same historical physical RHS, restricted to its downward, positive domain."""
    state = np.asarray(state, dtype=float)
    if state.shape != (3,) or not np.isfinite(state).all():
        raise ValueError("Expected finite [height, downward speed, mass]")
    if state[1] < 0 or state[2] <= 0:
        raise ValueError("The downward-entry model requires v >= 0 and m > 0")
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        value = np.asarray(particle_model(t, state), dtype=float)
    if not np.isfinite(value).all():
        raise FloatingPointError("Non-finite model derivative")
    return value


def invariant_residual(states):
    """Residual of the exact positive-domain mass/height/speed identity.

    A small residual is a numerical consistency diagnostic, not physical validation.
    """
    values = np.asarray(states, dtype=float)
    h, v, mass = values.T
    kappa = LAMBDA / (C_D * L_V)
    return np.log(mass / mass[0]) - kappa * (
        (v * v - v[0] * v[0]) / 2.0 + G * (h - h[0])
    )
