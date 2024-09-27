#!/usr/bin/env python
# coding: utf-8

import os
import sys
from glob import glob
from pathlib import Path

import fiona
import geopandas
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import snail
import snail.io
from shapely.geometry import LineString, Point, Polygon
from tqdm.auto import tqdm

# read in file that contains metadata (inc. column names for cost)
base_path = Path(sys.argv[1])
main_data_path = base_path / "processed_data"
meta = pd.read_csv(
    main_data_path / "networks/network_layers_hazard_intersections_details.csv"
)

# source: JEM/scripts/001-create-and-intersect-wgrid.py
"""
Usage
-----

To run on SoGE cluster:

    python compute-grid-value.py /soge-home/projects/mistral/jamaica-ccri

To run locally with a relative path, ensure all inputs are available, then:

    python compute-grid-value.py ./data/dir

"""


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


def create_grid(grid_id_tiff_path, boundary_file_path, empty_grid_path):
    bounds = gpd.read_file(boundary_file_path, layer="jamaica")

    # expand grid by buffer
    buffer = 100  # in meters (based on crs of bounds)
    boundary_box = bounds.bounds
    minx, miny, maxx, maxy = boundary_box.values[0]
    minx -= buffer
    miny -= buffer
    maxx += buffer
    maxy += buffer

    # cell side length in meters (based on crs of bounds)
    cell_length = 1000

    # determine grid bounding box to fit an integer number of grid cells in each dimension
    i, minx, maxx = harmonise_grid(minx, maxx, cell_length)
    j, miny, maxy = harmonise_grid(miny, maxy, cell_length)

    # create grid as TIFF and save to disk
    os.makedirs(os.path.dirname(empty_grid_path), exist_ok=True)
    command = f"gdal_create -outsize {i} {j} -a_srs EPSG:3448 -a_ullr {minx} {miny} {maxx} {maxy} {empty_grid_path}"  # note: projection is hard-coded
    os.system(command)

    # Load the .tiff file
    with rasterio.open(empty_grid_path) as src:
        # Extract grid data
        grid_ids = np.arange(src.width * src.height).reshape((src.height, src.width))
        grid_transform = src.transform
        grid_width = src.width
        grid_height = src.height

    # Export grid_ids as .tiff file
    with rasterio.open(
        grid_id_tiff_path,
        "w",
        driver="GTiff",
        height=grid_height,
        width=grid_width,
        count=1,
        dtype=rasterio.uint32,
        crs=src.crs,
        transform=grid_transform,
    ) as dst:
        dst.write(grid_ids.astype(rasterio.uint32), 1)


def read_grid(grid_id_tiff_path):
    hazard_paths = [grid_id_tiff_path]
    hazard_files = pd.DataFrame({"path": hazard_paths})
    hazard_files["key"] = [Path(path).stem for path in hazard_paths]
    hazard_files, grids = snail.io.extend_rasters_metadata(hazard_files)
    grid = grids[0]
    return grid, hazard_files


def intersect_grid_with_nodes(
    base_path, grid, grid_file, network_file, nodes_layer, output_file
):
    input_file = base_path / network_file

    points_df = gpd.read_file(input_file, layer=nodes_layer)
    points_df = points_df.to_crs("EPSG:3448")

    grid_intersections_p = snail.intersection.apply_indices(
        points_df, grid, index_i="i_0", index_j="j_0"
    )
    grid_intersections_p = snail.io.associate_raster_files(
        grid_intersections_p, grid_file
    )

    grid_intersections_p.rename(
        columns={"jamaica_1km_grid_ids": "grid_ids"}, inplace=True
    )
    grid_intersections_p.to_file(
        output_file,
        layer=nodes_layer,
        driver="gpkg",
    )


def intersect_grid_with_edges(
    base_path, grid, grid_file, network_file, edges_layer, output_file
):
    input_file = base_path / network_file

    edges = gpd.read_file(input_file, layer=edges_layer)
    edges = edges.to_crs("EPSG:3448")

    edges = snail.intersection.prepare_linestrings(edges)
    grid_intersections_l = snail.intersection.split_linestrings(edges, grid)
    grid_intersections_l["length_m"] = grid_intersections_l["geometry"].length
    grid_intersections_l = snail.intersection.apply_indices(
        grid_intersections_l, grid, index_i="i_0", index_j="j_0"
    )
    grid_intersections_l = snail.io.associate_raster_files(
        grid_intersections_l, grid_file
    )

    grid_intersections_l.rename(
        columns={"jamaica_1km_grid_ids": "grid_ids"}, inplace=True
    )
    grid_intersections_l.to_file(
        output_file,
        layer=edges_layer,
        driver="gpkg",
    )


def intersect_grid_with_polygons(
    base_path, grid, grid_file, network_file, polygon_layer, output_file
):
    input_file = base_path / network_file

    polygons_df = gpd.read_file(input_file, layer=polygon_layer)
    polygons_df = polygons_df.to_crs("EPSG:3448")

    polygons_df = snail.intersection.prepare_polygons(polygons_df)

    polygons_df = snail.intersection.split_polygons(polygons_df, grid)

    grid_intersections_poly = snail.intersection.apply_indices(
        polygons_df, grid, index_i="i_0", index_j="j_0"
    )
    grid_intersections_poly = snail.io.associate_raster_files(
        grid_intersections_poly, grid_file
    )

    grid_intersections_poly.rename(
        columns={"jamaica_1km_grid_ids": "grid_ids"}, inplace=True
    )
    grid_intersections_poly.to_file(
        output_file,
        layer=polygon_layer,
        driver="gpkg",
    )


# source: JEM/scripts/004-export-data-as-raster-grid.py


def setup_grid(tiff_file_path):
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
    with rasterio.open(output_path, "w", **output_kwargs) as dst:
        dst.write(output_grid, 1)


def read_csvs(csv_dir):
    for csv in tqdm(glob(str(csv_dir / "*.csv"))):
        yield pd.read_csv(csv)


def aggregate_to_grid(dfs, varname, shape, dtype):
    output_grid_flat = np.zeros(
        output_kwargs["height"] * output_kwargs["width"], dtype=dtype
    )

    for df in dfs:
        value = df[varname].sum()
        grid_id = df.loc[0, "grid_id"]
        output_grid_flat[grid_id] = value

    return output_grid_flat.reshape(shape)


def compute_mean_cost(df, min_cost, max_cost):
    df["cost_avg"] = df[[min_cost, max_cost]].mean(axis=1)
    return df


output_paths = []

ids_tiff = base_path / "results/grid_failures/jamaica_1km_grid_ids.tiff"
grid_ids, output_kwargs = setup_grid(ids_tiff)


base_path = Path(base_path)
# input paths
boundary_file_path = base_path / "processed_data/boundaries/jamaica.gpkg"
# output paths
empty_grid_path = base_path / "results/grid_failures/jamaica_1km_empty_grid.tiff"
grid_id_tiff_path = base_path / "results/grid_failures/jamaica_1km_grid_ids.tiff"
create_grid(grid_id_tiff_path, boundary_file_path, empty_grid_path)
grid, grid_file = read_grid(grid_id_tiff_path)


for _, row in meta.iterrows():
    data_path = Path(main_data_path / row["path"])
    base_name = data_path.stem
    output_path = Path(base_path / f"results/grid_failures/{base_name}_wgrid_ids.gpkg")
    layer = row["asset_layer"]
    data = gpd.read_file(data_path, layer=layer)

    if data.geometry.geom_type[0] == "Point":
        intersect_grid_with_nodes(
            base_path, grid, grid_file, data_path, layer, output_path
        )
    elif (
        data.geometry.geom_type[0] == "LineString"
        or data.geometry.geom_type[0] == "MultiLineString"
    ):
        intersect_grid_with_edges(
            base_path, grid, grid_file, data_path, layer, output_path
        )
    elif (
        data.geometry.geom_type[0] == "Polygon"
        or data.geometry.geom_type[0] == "MultiPolygon"
    ):
        intersect_grid_with_polygons(
            base_path, grid, grid_file, data_path, layer, output_path
        )
    else:
        print("unknown geometry type: ", data.geom_type[0], type(data.geom_type[0]))

    output_paths.append(output_path)


meta["output_path"] = output_paths

# meta.to_file(main_data_path / "networks/network_layers_hazard_intersections_details.csv")

for sector in meta["sector"].unique():
    sector_data = meta[meta["sector"] == sector]

    dtype = "float32"
    output_grid_flat = np.zeros(
        output_kwargs["height"] * output_kwargs["width"], dtype=dtype
    )

    for _, row in sector_data.iterrows():
        data_path = Path(main_data_path / row["output_path"])
        layer = row["asset_layer"]
        data = gpd.read_file(data_path, layer=layer)
        if (
            row["asset_mean_cost_column"] is not np.nan
        ):  # note: alternatively can apply this for water assets only
            cost_column = row["asset_mean_cost_column"]
        else:
            compute_mean_cost(
                data, row["asset_min_cost_column"], row["asset_max_cost_column"]
            )
            cost_column = "cost_avg"
        if layer == "nodes":
            for _, node in data.iterrows():
                output_grid_flat[node.grid_ids] += node[cost_column]
        elif layer == "edges":
            for _, edge in data.iterrows():
                output_grid_flat[edge.grid_ids] += edge.length_m * edge[cost_column]
        elif layer == "areas":
            for _, edge in data.iterrows():
                output_grid_flat[edge.grid_ids] += edge[cost_column]

    output_grid = output_grid_flat.reshape(grid_ids.shape)
    output_kwargs["dtype"] = dtype
    out_tiff = base_path / f"results/grid_failures/{sector}_avg_cost_grid.tiff"
    write_grid(out_tiff, output_grid, output_kwargs)
