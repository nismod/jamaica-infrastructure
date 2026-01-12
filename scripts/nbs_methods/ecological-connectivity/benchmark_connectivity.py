#!/usr/bin/env python3
"""
Simple benchmark harness for connectivity.landscape_connectivity.

This script intentionally keeps the data-loading flow close to Robynlibrary
so the timing reflects a realistic workload.
"""

from __future__ import annotations

import argparse
import os
import statistics
import time
from pathlib import Path

import numpy as np
import rasterio

import connectivity


def _load_condition_layer(path: Path, crop: tuple[int, int] | None) -> np.ndarray:
    with rasterio.open(path) as src:
        data = src.read(1)

    if crop:
        rows, cols = crop
        data = data[:rows, :cols]

    return data


def _build_land_mask(condition_layer: np.ndarray) -> np.ndarray:
    land_mask = np.zeros_like(condition_layer, dtype=bool)
    land_mask[condition_layer != 0] = True
    return land_mask


def run_benchmark(
    raster_path: Path,
    iterations: int,
    n_processes: int,
    lambda_parameter: float,
    gen_mode: str,
    number_of_generations: int,
    crop: tuple[int, int] | None,
) -> None:
    condition_layer = _load_condition_layer(raster_path, crop)
    land_array = _build_land_mask(condition_layer)

    print(f"Raster shape: {condition_layer.shape}, dtype={condition_layer.dtype}")
    print(
        f"Parameters: processes={n_processes}, lambda={lambda_parameter}, "
        f"mode={gen_mode}, generations={number_of_generations}"
    )
    if crop:
        print(f"Cropped to first {crop[0]} rows x {crop[1]} cols")

    timings: list[float] = []
    for iteration in range(1, iterations + 1):
        start = time.perf_counter()
        value = connectivity.landscape_connectivity(
            condition_layer,
            n_processes,
            land_array,
            lambda_parameter,
            gen_mode,
            number_of_generations,
            draw_plots=False,
        )
        elapsed = time.perf_counter() - start
        timings.append(elapsed)
        print(f"Iteration {iteration}: {elapsed:.2f}s (landscape connectivity={value:.4f})")

    if len(timings) == 1:
        summary = f"{timings[0]:.2f}s"
    else:
        summary = (
            f"min={min(timings):.2f}s median={statistics.median(timings):.2f}s "
            f"mean={statistics.mean(timings):.2f}s"
        )

    print(f"Benchmark complete: {summary}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark connectivity.py performance.")
    parser.add_argument(
        "--raster",
        type=Path,
        default=Path("Connectivity_inputs/landcover_condition_baseline_resampled.tif"),
        help="Condition raster to load for benchmarking.",
    )
    parser.add_argument("--iterations", type=int, default=1, help="Number of timed runs.")
    parser.add_argument(
        "--processes",
        type=int,
        default=os.cpu_count(),
        help="Number of worker processes to use.",
    )
    parser.add_argument(
        "--lambda-parameter",
        type=float,
        default=5.0,
        help="Lambda parameter passed through to landscape_connectivity.",
    )
    parser.add_argument(
        "--gen-mode",
        choices=("one_generation", "multi-generation"),
        default="one_generation",
        help="Generations mode used when computing connectivity.",
    )
    parser.add_argument(
        "--number-of-generations",
        type=int,
        default=1,
        help="Number of generations for the dispersal kernel.",
    )
    parser.add_argument(
        "--crop",
        type=str,
        help="Optional crop as ROWSxCOLS (e.g. 500x500) to benchmark a subset.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    crop = None
    if args.crop:
        try:
            row_str, col_str = args.crop.lower().split("x", 1)
            crop = (int(row_str), int(col_str))
        except ValueError as exc:
            raise SystemExit(f"Invalid crop '{args.crop}', expected ROWSxCOLS") from exc

    run_benchmark(
        raster_path=args.raster,
        iterations=args.iterations,
        n_processes=args.processes,
        lambda_parameter=args.lambda_parameter,
        gen_mode=args.gen_mode,
        number_of_generations=args.number_of_generations,
        crop=crop,
    )


if __name__ == "__main__":
    main()
