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
                                        "Flood Extent Rasters"
                                    )
    mangrove_output_path = os.path.join(
                                        incoming_data_path,
                                        "forces_of_nature_project_data",
                                        "Flood Extent Rasters",
                                        "merge_results"
                                    )
    if os.path.exists(mangrove_output_path) == False:
        os.mkdir(mangrove_output_path)

    return_periods = [25,100,500]
    mangrove_cases = ["NoMg_Ras","Mg_Ras"] 
    buffer_distance = 50
    flood_points = []
    for rp in return_periods:
        for mg in mangrove_cases:
            flood_file = f"RP_{rp}_{mg}"
            fpath_in = os.path.join(mangrove_input_path, f"{flood_file}.tif")
            fpath_out = os.path.join(mangrove_output_path, f"{flood_file}_reproject.tif")
            raster_rewrite(fpath_in,fpath_out)

            outCSVName = os.path.join(mangrove_output_path, f"{flood_file}_reproject.csv")
            subprocess.run(["gdal2xyz.py", '-csv', fpath_out, outCSVName])

            # Load points and convert to geodataframe with coordinates
            load_points = pd.read_csv(outCSVName, header=None, names=[
                                      "x", "y", flood_file], index_col=None)
            os.remove(outCSVName)
            os.remove(fpath_out)
            geometry = [Point(xy) for xy in zip(load_points.x, load_points.y)]
            load_points = gpd.GeoDataFrame(
                                            load_points, 
                                            crs=f"epsg:{epsg_jamaica}", 
                                            geometry=geometry
                                        )
            load_points["geometry"] = load_points.geometry.buffer(buffer_distance)
            if len(flood_points) > 0:
                load_points = load_points.drop(['x', 'y'], axis=1)
                matches = gpd.sjoin(flood_points,load_points, how="inner", op='intersects').reset_index()
                flood_points = pd.merge(flood_points,matches[['x','y',flood_file]],how="left",on=["x","y"])
            else:
                flood_points = load_points.copy()

    flood_points = flood_points.drop("geometry",axis=1)
    flood_points.to_csv(os.path.join(mangrove_output_path, "forces_to_nature_flood_results_raw.csv"),index=False)

    flood_points = pd.read_csv(os.path.join(mangrove_output_path, "forces_to_nature_flood_results_raw.csv"))
    flood_points["totals"] = 0
    for rp in return_periods:
        for mg in mangrove_cases:
            flood_points[f"RP_{rp}_{mg}"] = np.where(flood_points[f"RP_{rp}_{mg}"] < 0, 0, flood_points[f"RP_{rp}_{mg}"])
            flood_points["totals"] += flood_points[f"RP_{rp}_{mg}"]

    floods_filtered = flood_points[flood_points["totals"] > 0]
    flood_points.drop("totals",axis=1).to_csv(os.path.join(mangrove_output_path, "forces_to_nature_flood_results_masked.csv"),index=False)
    floods_filtered.drop("totals",axis=1).to_csv(os.path.join(mangrove_output_path, "forces_to_nature_flood_results_filtered.csv"),index=False)

if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
