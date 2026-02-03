"""
Create a raster grid encompassing a given boundary and write to disk (tiff) with
integer IDs as pixel values.
"""


import logging

import click
import geopandas as gpd
import numpy as np
import rasterio


def harmonise_grid(
    minimum: float, maximum: float, cell_length_meters: float
) -> tuple[int, float, float]:
    """
    Grow grid dimensions to encompass whole number of `cell_length_meters`

    Args:
        minimum: Minimum dimension value
        maximum: Maximum dimension value
        cell_length_meters: Length of cell side

    Returns:
        Number of cells
        Adjusted minimum
        Adjusted maximum
    """
    assert maximum > minimum

    span: float = maximum - minimum
    n_cells: int = int(np.ceil(span / cell_length_meters))
    delta: float = n_cells * cell_length_meters - span
    buffer: float = delta / 2

    return n_cells, minimum - buffer, maximum + buffer


@click.command()
@click.version_option("1.0")
@click.option(
    "--boundary-path", "-b", required=True, help="Path to country boundary file. Should be readable by geopandas.",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, readable=True),
)
@click.option(
    "--cell-length-meters", "-c", required=True, type=float,
    help="Width and height of grid spacing. In meters."
)
@click.option(
    "--boundary-buffer-meters", "-b", required=True, type=float,
    help="Approximate buffer to extend grid beyond boundary by. In meters.", 
)
@click.option(
    "--output-path", "-o", required=True, help="Path to write the output grid file",
    type=click.Path(exists=False, dir_okay=False, file_okay=True, readable=True),
)
def create_grid(
    boundary_path: str,
    cell_length_meters: float,
    boundary_buffer_meters: float,
    output_path: str,
) -> None:
    """
    Given a 2D boundary, generate a raster grid encompassing the provided
    boundary and write to disk as GeoTIFF.

    Args:
        boundary_path: Path to geospatial file containing single geometry. The
            output grid will be the bounding box of this geometry.
        cell_length_meters: Length of each output grid cell side.
        boundary_buffer_meters: Distance to buffer the boundary geometry by.
        output_path: Path to write the created grid to as GeoTIFF.
    """

    assert cell_length_meters > 0
    assert boundary_buffer_meters > 0

    boundary: gpd.GeoDataFrame = gpd.read_file(boundary_path)
    boundary = boundary.to_crs(epsg=3448)  # ensure we are working in meters
    bounding_geometry, = boundary.loc[:, "geometry"]

    minx, miny, maxx, maxy = bounding_geometry.bounds
    minx -= boundary_buffer_meters
    miny -= boundary_buffer_meters
    maxx += boundary_buffer_meters
    maxy += boundary_buffer_meters

    # determine grid bounding box to fit an integer number of grid cells in each dimension
    i, minx, maxx = harmonise_grid(minx, maxx, cell_length_meters)
    j, miny, maxy = harmonise_grid(miny, maxy, cell_length_meters)

    transform = rasterio.Affine(cell_length_meters, 0, minx, 0, cell_length_meters, miny)

    logging.info(f"Writing grid to disk with {cell_length_meters=}, {i=} and {j=}")
    logging.info(f"Transform:\n{transform}")

    # export grid_ids as .tiff file
    with rasterio.open(
        output_path,
        "w",
        driver="GTiff",
        height=j,
        width=i,
        count=1,
        dtype="int32",
        crs=boundary.crs,
        transform=transform
    ) as dataset:
        dataset.write(np.zeros((j, i)), 1)


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    create_grid()
