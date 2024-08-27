"""Create the pipelines asset data for Jamaica
    Add the final asset data into a geopackage
"""

import os

import pandas as pd
import geopandas as gpd
import numpy as np
import zipfile
import io
from fiona.io import ZipMemoryFile
from fiona import BytesCollection
from preprocess_utils import *
from tqdm import tqdm

tqdm.pandas()


def gdf_geom_clip(gdf_in, clip_geom):
    """Filter a dataframe to contain only features within a clipping geometry

    Parameters
    ---------
    gdf_in
        geopandas dataframe to be clipped in
    clip_geom
        shapely geometry of that we need to intersect

    Returns
    -------
    filtered dataframe
    """
    return gdf_in.loc[
        gdf_in["geometry"].apply(lambda x: x.intersects(clip_geom))
    ].reset_index(drop=True)


def main(config):
    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]
    epsg_jamaica = 3448

    hydrobasins_outputs = os.path.join(
        processed_data_path, "hydrobasins", "hydrobasins.gpkg"
    )

    jamaica_boundary = gpd.read_file(
        os.path.join(processed_data_path, "boundaries", "admin_boundaries.gpkg"),
        layer="admin0",
    )
    jamaica_boundary = jamaica_boundary.to_crs(epsg=epsg_jamaica)

    root_dir = os.path.join(incoming_data_path, "hydrobasins", "na")

    for root, dirs, files in os.walk(root_dir):
        for file in files:
            if file.endswith(".zip"):
                root_file = os.path.join(root, file)
                zipshp = io.BytesIO(open(root_file, "rb").read())
                with BytesCollection(zipshp.read()) as src:
                    crs = src.crs
                    basins = gpd.GeoDataFrame.from_features(src, crs=crs)
                    basins = basins.to_crs(epsg=epsg_jamaica)
                    basins = gdf_geom_clip(basins, jamaica_boundary.geometry.values[0])
                    basins.to_file(
                        hydrobasins_outputs,
                        layer=file.replace(".zip", ""),
                        driver="GPKG",
                    )

                print(f"* Done with file {file}")


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
