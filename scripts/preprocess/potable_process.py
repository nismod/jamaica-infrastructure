"""Create the potable water asset data for Jamaica
    Add the final asset data into a geopackage
"""
import os

import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
from preprocess_utils import *


def main(config):
    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]
    epsg_jamaica = 3448
    potable_data_path = os.path.join(incoming_data_path, "water", "potable", "ss")

    # read in raw files

    # NOTE assumptions:
    # - pipes and storage tanks are assumed to not be susceptible to flooding/hurricane hazards and are therefore not provided with a cost or damage function
    # - water treatment plants are assigned a cost based on capacity: limited capacity data was supplemented by matching with system attributes from Parish water supply plans where available
    # - hurricane induced damage for pumping stations, intakes and wells are assumed to be the same as that for water treatment plants
    # - flood induced damage for intakes is assumed to be the same as that for pumping stations

    potable_facilities_NWC = pd.read_csv(
        os.path.join(potable_data_path, "potable_facilities_NWC_test.csv")
    )
    asset_population_served = pd.read_csv(
        os.path.join(potable_data_path, "final_3.csv")
    )
    potable_facilities_NWC = pd.merge(
        potable_facilities_NWC,
        asset_population_served[["OBJECTID", "path"]],
        on="OBJECTID",
    ).rename(columns={"Type": "asset_type"})

    # provide id
    potable_facilities_NWC["node_id"] = potable_facilities_NWC \
        .apply(lambda node: f"{node.asset_type}_{node.OBJECTID}", axis=1)

    # export as gpkg
    geometry = [
        Point(xy)
        for xy in zip(
            potable_facilities_NWC["lon"].round(4).values,
            potable_facilities_NWC["lat"].round(4).values,
        )
    ]
    potable_facilities_NWC = gpd.GeoDataFrame(
        potable_facilities_NWC, crs=f"EPSG:{epsg_jamaica}", geometry=geometry
    )
    potable_facilities_NWC.to_file(
        os.path.join(
            processed_data_path, "networks", "water", "potable_facilities_NWC.gpkg"
        ),
        layer="nodes",
        driver="GPKG",
    )


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
