"""
Create a raster grid encompassing a given boundary and write to disk (tiff) with
integer IDs as pixel values.
"""


from pathlib import Path
import sys

import geopandas as gpd
import numpy as np
import rasterio
import rasterio.crs


def harmonise_grid(
    minimum: float, maximum: float, cell_length: float
) -> tuple[int, float, float]:
    """
    Grow grid dimensions to encompass whole number of `cell_length`

    Args:
        minimum: Minimum dimension value
        maximum: Maximum dimension value
        cell_length: Length of cell side

    Returns:
        Number of cells
        Adjusted minimum
        Adjusted maximum
    """
    assert maximum > minimum

    span: float = maximum - minimum
    n_cells: int = int(np.ceil(span / cell_length))
    delta: float = n_cells * cell_length - span
    buffer: float = delta / 2

    return n_cells, minimum - buffer, maximum + buffer


def create_grid(
    boundary_path: Path,
    cell_length: float,
    boundary_buffer: float,
    output_grid: Path,
) -> None:
    """
    Given a 2D boundary, generate a raster grid encompassing the provided
    boundary and write to disk as GeoTIFF.

    Args:
        boundary_path: Path to geospatial file containing single geometry. The
            output grid will be the bounding box of this geometry.
        cell_length: Length of each output grid cell side. In units of boundary CRS.
        boundary_buffer: Distance to buffer the boundary geometry by. In units of boundary CRS.
        output_grid: Path to write the created grid to as GeoTIFF.
    """

    assert cell_length > 0
    assert boundary_buffer > 0

    boundary: gpd.GeoDataFrame = gpd.read_file(boundary_path)
    bounding_geometry, = boundary.loc[:, "geometry"]

    minx, miny, maxx, maxy = bounding_geometry.bounds
    minx -= boundary_buffer
    miny -= boundary_buffer
    maxx += boundary_buffer
    maxy += boundary_buffer

    # determine grid bounding box to fit an integer number of grid cells in each dimension
    i, minx, maxx = harmonise_grid(minx, maxx, cell_length)
    j, miny, maxy = harmonise_grid(miny, maxy, cell_length)

    # export grid_ids as .tiff file
    with rasterio.open(
        output_grid,
        "w",
        driver="GTiff",
        height=j,
        width=i,
        count=1,
        dtype="int32",
        crs=boundary.crs,
        transform=rasterio.Affine(cell_length, 0, minx, 0, cell_length, miny)
    ) as dataset:
        dataset.write(np.zeros((j, i)), 1)


if __name__ == "__main__":
    try:
        boundary_path = Path(sys.argv[1])
        cell_length = float(sys.argv[2])
        buffer = float(sys.argv[3])
        output_grid = Path(sys.argv[4])
    except IndexError:
        print(
            "Usage: \npython generate_grid.py <input_boundary_path> <cell_length> <boundary_buffer> <output_grid_path>\n"
            "Note that <cell_length> and <boundary_buffer> must be in units of the boundary CRS."
        )
        sys.exit()

    create_grid(boundary_path, cell_length, buffer, output_grid)
