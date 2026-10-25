# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Fengrui Zhang, Yiye Huang, Zening Wang
"""Run bounded demonstrations or explicitly sized parameter scans."""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
import matplotlib

from space_dust.integrate import integrate
from space_dust.model import initial_state
from space_dust.plotting import plot_scan, plot_trajectory


def simulate(args, radius, speed):
    result = integrate(initial_state(radius, speed, args.altitude_km),
                       t_max=args.t_max, dt=args.dt,
                       mass_fraction_cutoff=args.mass_cutoff,
                       max_steps=args.max_steps)
    row = result.summary()
    row.update(initial_radius_mm=radius, initial_speed_kms=speed,
               initial_altitude_km=args.altitude_km, requested_dt_s=args.dt,
               requested_t_max_s=args.t_max)
    return result, row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for name in ["demo", "scan"]:
        sub = subparsers.add_parser(name)
        sub.add_argument("--output", type=Path, default=Path("results") / name)
        sub.add_argument("--dt", type=float, default=0.01)
        sub.add_argument("--t-max", type=float, default=100.0)
        sub.add_argument("--altitude-km", type=float, default=100.0)
        sub.add_argument("--mass-cutoff", type=float, default=1e-12,
                         help="Relative numerical mass threshold; not proof of complete ablation")
        sub.add_argument("--max-steps", type=int, default=200_000)
        if name == "demo":
            sub.add_argument("--radius-mm", type=float, default=0.1)
            sub.add_argument("--speed-kms", type=float, default=12.0)
        else:
            sub.add_argument("--radii-mm", type=float, nargs="+", default=[0.05, 0.1])
            sub.add_argument("--speeds-kms", type=float, nargs="+", default=[12.0, 20.0])
            sub.add_argument("--allow-large-scan", action="store_true",
                             help="Explicitly allow more than 16 parameter combinations")
    args = parser.parse_args()
    if args.command == "scan":
        if len(set(args.radii_mm)) != len(args.radii_mm) or len(set(args.speeds_kms)) != len(args.speeds_kms):
            parser.error("Scan axes must not contain duplicates")
        if len(args.radii_mm) * len(args.speeds_kms) > 16 and not args.allow_large_scan:
            parser.error("More than 16 cases requires --allow-large-scan")
    if args.output.exists() and (not args.output.is_dir() or any(args.output.iterdir())):
        parser.error("Output must be new or empty; existing results are never overwritten")
    args.output.mkdir(parents=True, exist_ok=True)
    records = []
    if args.command == "demo":
        result, row = simulate(args, args.radius_mm, args.speed_kms)
        records.append(row)
        np.savez_compressed(args.output / "trajectory.npz", time_s=result.time_s, states=result.states)
        plot_trajectory(result, args.output / "trajectory.png")
    else:
        for speed in args.speeds_kms:
            for radius in args.radii_mm:
                _, row = simulate(args, radius, speed)
                records.append(row)
        with (args.output / "scan.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
        plot_scan(records, args.radii_mm, args.speeds_kms, args.output / "scan.png")
    report = {
        "solver": "guarded_rk4; original physical RHS; not the historical notebook execution",
        "environment": {"python": sys.version.split()[0], "numpy": np.__version__,
                        "matplotlib": matplotlib.__version__},
        "mass_fraction_cutoff": args.mass_cutoff,
        "records": records,
        "warning": "time_limit is not ground arrival; remaining mass is not a ground survival rate",
    }
    (args.output / "summary.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    counts = {status: sum(r["termination"] == status for r in records)
              for status in sorted({r["termination"] for r in records})}
    print(json.dumps({"cases": len(records), "termination_counts": counts}))
    return 2 if any(r["termination"] in {"numerical_failure", "step_limit"} for r in records) else 0


if __name__ == "__main__":
    raise SystemExit(main())
