# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Fengrui Zhang, Yiye Huang, Zening Wang
"""Domain-guarded RK4 with explicit endpoint classifications.

Step halving rejects invalid stages; it is not an error-controlled adaptive solver.
Ground/cutoff events are linearly interpolated within a valid accepted step.
"""
from dataclasses import dataclass
import numpy as np

from .model import invariant_residual, positive_rhs, radius_from_mass


@dataclass
class Trajectory:
    time_s: np.ndarray
    states: np.ndarray
    status: str
    detail: str
    accepted_steps: int
    rejected_steps: int

    def summary(self):
        start, end = self.states[0], self.states[-1]
        fraction = float(end[2] / start[2])
        return {
            "termination": self.status,
            "termination_detail": self.detail,
            "time_s": float(self.time_s[-1]),
            "height_km": float(end[0] / 1e3),
            "downward_speed_m_s": float(end[1]),
            "mass_kg": float(end[2]),
            "mass_mg": float(end[2] * 1e6),
            "radius_mm": float(radius_from_mass(end[2]) * 1e3),
            "end_mass_percent": 100.0 * fraction,
            "end_radius_percent": 100.0 * fraction ** (1.0 / 3.0),
            "ground_mass_percent": 100.0 * fraction if self.status == "ground_reached" else None,
            "accepted_steps": self.accepted_steps,
            "rejected_steps": self.rejected_steps,
            "maximum_abs_invariant_residual": float(np.max(np.abs(invariant_residual(self.states)))),
            "interpretation": (
                "End-of-run remaining mass is not ground survival unless termination is ground_reached. "
                "mass_cutoff is a numerical threshold, not demonstrated complete ablation."
            ),
        }


def _evaluate(rhs, time, state):
    if not np.isfinite(state).all() or state[1] < 0 or state[2] <= 0:
        raise ValueError("Invalid intermediate RK4 stage")
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        value = np.asarray(rhs(time, state), dtype=float)
    if value.shape != (3,) or not np.isfinite(value).all():
        raise FloatingPointError("Invalid derivative")
    return value


def _rk4_step(rhs, time, state, step):
    k1 = _evaluate(rhs, time, state)
    k2 = _evaluate(rhs, time + step / 2, state + step * k1 / 2)
    k3 = _evaluate(rhs, time + step / 2, state + step * k2 / 2)
    k4 = _evaluate(rhs, time + step, state + step * k3)
    candidate = state + step * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    if not np.isfinite(candidate).all() or candidate[1] < 0 or candidate[2] <= 0:
        raise ValueError("Invalid completed RK4 state")
    if candidate[2] > state[2] * (1 + 1e-12):
        raise ValueError("Mass increase is outside the ablation model")
    return candidate


def integrate(initial, t_max=100.0, dt=0.01, *, rhs=positive_rhs,
              mass_fraction_cutoff=1e-12, max_steps=200_000,
              max_halvings=24, min_dt=1e-10):
    state = np.asarray(initial, dtype=float).copy()
    if state.shape != (3,) or not np.isfinite(state).all():
        raise ValueError("Initial state must be finite [h, v, m]")
    if state[0] < 0 or state[1] < 0 or state[2] <= 0:
        raise ValueError("Initial state requires h >= 0, v >= 0, m > 0")
    if not all(np.isfinite(x) and x > 0 for x in [t_max, dt, min_dt]):
        raise ValueError("t_max, dt, min_dt must be positive and finite")
    if not 0 < mass_fraction_cutoff < 1 or max_steps < 1 or max_halvings < 0:
        raise ValueError("Invalid cutoff or iteration limits")
    cutoff = state[2] * mass_fraction_cutoff
    times, states = [0.0], [state.copy()]
    accepted = rejected = 0
    status, detail = "time_limit", "Requested time horizon reached without a ground/cutoff event"
    current_step = dt
    if state[0] == 0:
        status, detail = "ground_reached", "Initial state is at the ground"
    else:
        while times[-1] < t_max:
            if accepted >= max_steps:
                status, detail = "step_limit", "Maximum accepted-step budget reached"
                break
            time = times[-1]
            trial = min(current_step, t_max - time)
            candidate = None
            error = ""
            for _ in range(max_halvings + 1):
                try:
                    candidate = _rk4_step(rhs, time, state, trial)
                    break
                except (ArithmeticError, ValueError) as exc:
                    error = str(exc)
                    rejected += 1
                    trial /= 2
                    if trial < min_dt:
                        break
            if candidate is None:
                status, detail = "numerical_failure", "No valid RK4 stage sequence: " + error
                break
            events = []
            if candidate[0] <= 0:
                events.append((state[0] / (state[0] - candidate[0]), "ground_reached"))
            if candidate[2] <= cutoff:
                events.append(((state[2] - cutoff) / (state[2] - candidate[2]), "mass_cutoff"))
            accepted += 1
            if events:
                fraction, status = min(events, key=lambda item: item[0])
                state = state + fraction * (candidate - state)
                state[0 if status == "ground_reached" else 2] = 0.0 if status == "ground_reached" else cutoff
                times.append(time + fraction * trial)
                states.append(state.copy())
                detail = "Event linearly interpolated within a valid accepted step"
                break
            state = candidate
            times.append(min(time + trial, t_max))
            states.append(state.copy())
            current_step = min(dt, 2 * trial)
    return Trajectory(np.asarray(times), np.asarray(states), status, detail, accepted, rejected)
