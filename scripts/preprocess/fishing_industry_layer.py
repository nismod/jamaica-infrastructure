"""Process a fishery layer to assign farm sizes to farm areas
    Gap fill the missing farm size information 
    By fitting a linear curve between known farms sizes and given areas
"""

import sys
import os

import geopandas as gpd
import pandas as pd
import numpy as np
from preprocess_utils import *
from tqdm import tqdm

tqdm.pandas()


def main(config):
    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]
    epsg_jamaica = 3448

    fishing_locations = gpd.read_file(
        os.path.join(
            incoming_data_path, "nsdmb", "GWP_Jamaica_NSP_Master_Geodatabase_v01.gdb"
        ),
        layer="aqua_farms",
    )
    fishing_locations = fishing_locations.to_crs(epsg=epsg_jamaica)
    prod_locations = fishing_locations[fishing_locations["Size_Farm"] > 0]
    polyfit = np.polyfit(
        prod_locations["Shape_Area"], prod_locations["Size_Farm"], deg=1
    )
    noprod_locations = fishing_locations[fishing_locations["Size_Farm"] == 0]
    noprod_locations["Size_Farm"] = (
        polyfit[0] * noprod_locations["Shape_Area"] + polyfit[1]
    )

    fishing_locations = gpd.GeoDataFrame(
        pd.concat([prod_locations, noprod_locations], axis=0, ignore_index=True),
        geometry="geometry",
        crs=f"EPSG:{epsg_jamaica}",
    )
    fishing_locations["farm_id"] = fishing_locations.index.values.tolist()
    fishing_locations.to_file(
        os.path.join(processed_data_path, "land_type_and_use", "aqua_farms.gpkg"),
        layer="areas",
        driver="GPKG",
    )


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
