# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Fengrui Zhang, Yiye Huang, Zening Wang
"""Plots with endpoint labels that distinguish remaining mass from ground arrival."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .model import radius_from_mass


def plot_trajectory(result, destination):
    states, times = result.states, result.time_s
    mass_fraction = states[:, 2] / states[0, 2]
    fig, axes = plt.subplots(4, 1, figsize=(8, 10), sharex=True)
    values = [states[:, 0] / 1e3, states[:, 1] / 1e3,
              mass_fraction * 100, np.cbrt(mass_fraction) * 100]
    labels = ["Altitude (km)", "Downward speed (km/s)", "Mass remaining (%)", "Radius remaining (%)"]
    for ax, value, label in zip(axes, values, labels):
        ax.plot(times, value, color="#315c76")
        ax.set_ylabel(label)
        ax.grid(alpha=0.25)
    axes[0].set_title("Guarded RK4 demonstration | termination: " + result.status)
    axes[-1].set_xlabel("Time (s)")
    fig.tight_layout()
    fig.savefig(destination, dpi=160)
    plt.close(fig)


def plot_scan(records, radii, speeds, destination):
    from matplotlib.colors import BoundaryNorm, ListedColormap
    statuses = ["ground_reached", "mass_cutoff", "time_limit", "step_limit", "numerical_failure"]
    mass = np.array([r["end_mass_percent"] for r in records]).reshape(len(speeds), len(radii))
    flags = np.array([statuses.index(r["termination"]) for r in records]).reshape(mass.shape)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.3), constrained_layout=True)
    image = axes[0].imshow(mass, origin="lower", aspect="auto", vmin=0, vmax=100, cmap="viridis")
    fig.colorbar(image, ax=axes[0], label="Mass remaining at termination (%)")
    cmap = ListedColormap(["#4c956c", "#e9c46a", "#8da0cb", "#b3b3b3", "#e76f51"])
    status_image = axes[1].imshow(flags, origin="lower", aspect="auto", cmap=cmap,
                                  norm=BoundaryNorm(np.arange(-0.5, 5.5), 5))
    bar = fig.colorbar(status_image, ax=axes[1], ticks=range(5))
    bar.ax.set_yticklabels(statuses)
    for ax in axes:
        ax.set_xticks(range(len(radii)), [f"{x:g}" for x in radii])
        ax.set_yticks(range(len(speeds)), [f"{x:g}" for x in speeds])
        ax.set_xlabel("Initial radius (mm)")
        ax.set_ylabel("Initial speed (km/s)")
    axes[0].set_title("End-of-run mass, not ground survival")
    axes[1].set_title("Termination classification")
    fig.savefig(destination, dpi=160)
    plt.close(fig)
