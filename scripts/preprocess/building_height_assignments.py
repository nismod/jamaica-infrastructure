""" Assign heights to buildings in Jamaica
"""

import os
import glob
import geopandas as gpd
import pandas as pd
import numpy as np
import subprocess
from preprocess_utils import *
from tqdm import tqdm
tqdm.pandas()

def main(config):
    # Set global paths
    processed_data_path = config['paths']['data']
    results_data_path = config['paths']['results']

    # Create a path to store intersection outputs
    output_path = os.path.join(results_data_path,
                            "buildings_heights")
    if os.path.exists(output_path) == False:
        os.mkdir(output_path) 

    vector_details_csv = os.path.join(processed_data_path,"buildings","building_layers.csv")
    raster_details_csv = os.path.join(processed_data_path,"buildings","building_rasters.csv")

    run_intersections = True  # Set to True is you want to run this process
    if run_intersections is True:
        args = [
                "python",
                "vector_raster_intersections.py",
                f"{vector_details_csv}",
                f"{raster_details_csv}",
                f"{output_path}"
                ]
        print ("* Start the processing of Roads voronoi and population raster intersections")
        print (args)
        subprocess.run(args)

    print ("* Done with the processing of Roads voronoi and population raster intersections")

    """Post-processing the buildings-raster intersection results
    """
    building_id_column = "osm_id" # Building ID column
    raster_columns = pd.read_csv(raster_details_csv)["key"].values.tolist()
    # Read in intersection geoparquet
    intersections_file = [file for file in os.listdir(output_path) if file.endswith(".geoparquet")][0]
    if intersections_file:
        building_intersections = gpd.read_parquet(intersections_file)
        all_buildings = gpd.read_file(os.path.join(
                                processed_data_path,
                                "buildings",
                                "buildings_assigned_economic_activity.gpkg"),
                                layer="areas")
        for col in raster_columns:
            building_values = building_intersections[building_intersections[col] > 0]
            building_values = building_values.groupby(building_id_column)[col].mean().reset_index()
            all_buildings = pd.merge(all_buildings,building_values,how="left",on=[building_id_column]).fillna(0)
            del building_values
    
        del building_intersections

    all_buildings.to_file(
                        os.path.join(
                                processed_data_path,
                                "buildings",
                                "buildings_heights.gpkg"
                            ),
                        layer="areas",driver="GPKG")
    print("* Done with assigning heights to buildings")



if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
