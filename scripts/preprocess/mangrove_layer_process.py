"""Create a combined land use layer from TNC and Forestry land use layers
    Also combine a global mining areas land use layer 
"""
import sys
import os
import subprocess

import rasterio
import rioxarray
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import numpy as np
from preprocess_utils import *
from tqdm import tqdm
tqdm.pandas()

epsg_jamaica = 3448

def raster_rewrite(in_raster,out_raster):
    """Rewrite a raster to reproject and change no data value

    Parameters
        - in_raster - String name of input GeoTff file path
        - out_raster - String name of output GeoTff file path
        - nodata - Float value of data that is treated as no data

    Outputs
        Reproject and replace raster with nodata = -1
    """
    ds = rioxarray.open_rasterio(in_raster, mask_and_scale=True)
    ds = ds.rio.reproject(f"EPSG:{epsg_jamaica}")
    ds.rio.to_raster(out_raster)

    # os.remove(in_raster)
    # os.rename(out_raster,in_raster)

def main(config):
    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]

    mangrove_input_path = os.path.join(
                                        incoming_data_path,
                                        "forces_of_nature_project_data",
                                        "mangroves/commondata/nsdmd",
                                    )


    mangroves_file = gpd.read_file(os.path.join(mangrove_input_path,"mangroves.shp"))
    print (mangroves_file)
    mangroves_file = mangroves_file.to_crs(epsg=epsg_jamaica)
    mangroves_file.to_file(os.path.join(processed_data_path,"nbs","nsmdb-mangroves.gpkg"))

if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
