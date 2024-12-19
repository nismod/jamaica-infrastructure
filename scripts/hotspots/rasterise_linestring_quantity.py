"""
Given a quantity defined on a set of linestrings, this script rasterises the quantity on a grid.
"""

import logging
import sys

from affine import Affine
import geopandas as gpd
import numpy as np
import pandas as pd
import rioxarray as rio
import xarray as xr

import snail
import snail.intersection


if __name__ == "__main__":

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)

    try:
        edges_path = sys.argv[1]
        grid_path = sys.argv[2]
        output_path = sys.argv[3]
        variable_name = sys.argv[4]
        aggregation = sys.argv[5]
    except IndexError:
        print(
            "Usage:\n"
            f"python {sys.argv[0]} <edges> <grid> <output> <edges_variable_name> <aggregation>\n"
            "e.g.\n"
            f"python {sys.argv[0]} edges.gpkg grid.tif output.tif rerouting_loss sum"
        )
        sys.exit(1)

    logging.info("Read vector data")
    raw_edges = gpd.read_file(edges_path)
    assert "_i_index" not in raw_edges.columns
    assert "_j_index" not in raw_edges.columns

    edges: gpd.GeoDataFrame = raw_edges[~raw_edges.geometry.isna()]
    edges: gpd.GeoDataFrame = snail.intersection.prepare_linestrings(edges)
    assert {"LineString"} == set(edges.geometry.type)

    logging.info("Read raster grid")
    grid_metadata = snail.intersection.GridDefinition.from_raster(grid_path)
    grid: xr.DataArray = rio.open_rasterio(grid_path)

    logging.info(f"Splitting linestrings by grid, rasterising {variable_name}")
    splits: gpd.GeoDataFrame = snail.intersection.split_linestrings(edges, grid_metadata)
    splits: gpd.GeoDataFrame = snail.intersection.apply_indices(splits, grid_metadata, index_i="_i_index", index_j="_j_index")
    logging.info(f"Aggregating {variable_name} by {aggregation}")
    accumulated: pd.DataFrame = (
        splits.loc[:, [variable_name, "_i_index", "_j_index"]]
        .groupby(["_i_index", "_j_index"]).aggregate(aggregation)
    )
    arr: np.ndarray = np.zeros((grid_metadata.height, grid_metadata.width))
    arr[accumulated.index.get_level_values("_j_index"), accumulated.index.get_level_values("_i_index")] \
        = accumulated.loc[:, variable_name].values

    logging.info("Writing raster to disk")
    ds: xr.DataArray = (
        xr.DataArray(arr, name=variable_name, coords={"y": grid.y, "x": grid.x})
        .astype("float32")
        .rio.write_nodata(-1.0)
        .rio.write_transform(Affine(*grid_metadata.transform))
        .rio.write_crs(grid_metadata.crs.to_epsg())
    )
    ds.rio.to_raster(output_path)
