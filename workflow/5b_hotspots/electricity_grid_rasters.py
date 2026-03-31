"""
Aggregate grid cell disruption results and write to raster (GeoTIFFs).

Takes a directory of grid cell CSV files and aggregates specified variables
to the hotspots grid, outputting as GeoTIFF rasters.
"""

import logging

import click
import numpy as np
import pandas as pd
import rasterio
from glob import glob
from pathlib import Path
from tqdm.auto import tqdm


def setup_grid(tiff_file_path):
    """Read grid metadata for output raster creation."""
    with rasterio.open(tiff_file_path) as src:
        grid_ids = src.read(1)
        grid_transform = src.transform
        grid_width = src.width
        grid_height = src.height
        grid_crs = src.crs

    output_kwargs = {
        "driver": "GTiff",
        "count": 1,
        "width": grid_width,
        "height": grid_height,
        "crs": grid_crs,
        "transform": grid_transform,
        "compress": "lzw",
    }

    return grid_ids, output_kwargs


def write_grid(output_path, output_grid, output_kwargs):
    """Write a grid array to a GeoTIFF file."""
    logging.info(f"Writing raster to {output_path}")
    with rasterio.open(output_path, "w", **output_kwargs) as dst:
        dst.write(output_grid, 1)


def read_csvs(csv_dir, pattern):
    """Generator to read CSV files matching pattern."""
    csv_files = list(glob(str(csv_dir / pattern)))
    logging.info(f"Reading {len(csv_files)} CSV files")
    for csv in tqdm(csv_files, desc="Reading CSVs"):
        yield pd.read_csv(csv)


def aggregate_to_grid(dfs, varname, shape, dtype):
    """Aggregate a variable from dataframes to grid."""
    logging.info(f"Aggregating {varname} to grid")
    height, width = shape
    output_grid = np.zeros((height, width), dtype=dtype)

    for df in dfs:
        if df.empty or varname not in df.columns:
            continue
        # Group by cell coordinates and sum values
        grouped = df.groupby(["cell_index_y", "cell_index_x"])[varname].sum()
        for (y, x), value in grouped.items():
            if 0 <= y < height and 0 <= x < width:
                output_grid[y, x] = value

    return output_grid


@click.command()
@click.version_option("1.0.0")
@click.option(
    "--grid-path",
    "-g",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Path to the hotspots grid TIFF file.",
)
@click.option(
    "--disruption-dir",
    "-d",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True),
    help="Directory containing disruption CSV files.",
)
@click.option(
    "--loss-dir",
    "-l",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True),
    help="Directory containing loss CSV files.",
)
@click.option(
    "--output-population-path",
    "-p",
    required=True,
    type=click.Path(dir_okay=False, writable=True),
    help="Output path for population affected raster.",
)
@click.option(
    "--output-demand-path",
    "-m",
    required=True,
    type=click.Path(dir_okay=False, writable=True),
    help="Output path for demand affected raster.",
)
@click.option(
    "--output-gdp-path",
    "-o",
    required=True,
    type=click.Path(dir_okay=False, writable=True),
    help="Output path for GDP loss raster.",
)
def main(
    grid_path,
    disruption_dir,
    loss_dir,
    output_population_path,
    output_demand_path,
    output_gdp_path,
):
    """Aggregate disruption and loss data to raster grids."""
    logging.info(f"Reading grid from {grid_path}")
    grid_ids, output_kwargs = setup_grid(grid_path)

    disruption_dir = Path(disruption_dir)
    loss_dir = Path(loss_dir)

    # 1. Aggregate population affected
    logging.info("Processing population_affected")
    varname = "population_affected"
    dtype = "float32"
    dfs = read_csvs(disruption_dir, "disruption_*.csv")
    output_grid = aggregate_to_grid(dfs, varname, grid_ids.shape, dtype)
    output_kwargs["dtype"] = dtype
    write_grid(output_population_path, output_grid, output_kwargs)

    # 2. Aggregate demand affected
    logging.info("Processing demand_affected")
    varname = "demand_affected"
    dtype = "int32"
    dfs = read_csvs(disruption_dir, "disruption_*.csv")
    output_grid = aggregate_to_grid(dfs, varname, grid_ids.shape, dtype)
    output_kwargs["dtype"] = dtype
    write_grid(output_demand_path, output_grid, output_kwargs)

    # 3. Aggregate GDP loss
    logging.info("Processing loss_gdp")
    varname = "loss_gdp"
    dtype = "float32"
    dfs = read_csvs(loss_dir, "loss_*.csv")
    output_grid = aggregate_to_grid(dfs, varname, grid_ids.shape, dtype)
    output_kwargs["dtype"] = dtype
    write_grid(output_gdp_path, output_grid, output_kwargs)


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    main()
    logging.info("Complete")
