"""Create the irrigation asset data for Jamaica
    Add the final asset data into a geopackage
"""
import os

import pandas as pd
import geopandas as gpd
import numpy as np
from preprocess_utils import *


def main(config):
    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]
    epsg_jamaica = 3448
    irrigation_data_path = os.path.join(
        incoming_data_path, "water", "irrigation", "raw"
    )

    # Read in raw files, assign asset name for matching with cost/damage data, concat

    ## NOTE assumptions:
    # - reservoirs and dams (surface water sources) are not susceptible to flooding or hurricanes
    # - the construction cost of potable water well is the same as irrigation well
    # - the flood damage function for canals is the same as for irrigation pipelines


    ## currently ignoring reservoirs and dams
    # reservoirs = gpd.read_file(
    #     os.path.join(irrigation_data_path, "WATER SOURCES", "RESERVOIRS.shp")
    # )
    # dams = gpd.read_file(os.path.join(irrigation_data_path, "WATER SOURCES", "DAM.shp"))
    # micro_dams = gpd.read_file(
    #     os.path.join(irrigation_data_path, "WATER SOURCES", "MICRO_DAMS.shp")
    # )
    wells = gpd.read_file(
        os.path.join(irrigation_data_path, "Well Sites", "NIC_WELL_SITES_2020.shp")
    )
    canals = gpd.read_file(os.path.join(irrigation_data_path, "canal_network_NIC.shp"))[
        ["OBJECTID", "geometry"]
    ]
    pipelines = gpd.read_file(
        os.path.join(irrigation_data_path, "pipelines_network_NIC.shp")
    )[["OBJECTID", "geometry"]]

    wells["OBJECTID"] = np.linspace(0, wells.shape[0], wells.shape[0])
    wells = wells[["OBJECTID", "geometry"]]
    wells["asset_type"] = "well"
    wells["asset_type_cost_data"] = "well"
    wells["asset_type_flood_damage"] = "well"
    wells["asset_type_hurricane_damage"] = "na"

    canals["asset_type"] = "canal"
    canals["asset_type_cost_data"] = "irrigation_canal"
    canals["asset_type_flood_damage"] = "irrigation_canal"
    canals["asset_type_hurricane_damage"] = "irrigation_canal"

    pipelines["asset_type"] = "pipeline"
    pipelines["asset_type_cost_data"] = "irrigation_pipeline"
    pipelines["asset_type_flood_damage"] = "irrigation_canal"
    pipelines["asset_type_hurricane_damage"] = "na"

    nodes = gpd.GeoDataFrame(wells, crs=f"EPSG:{epsg_jamaica}")
    edges = gpd.GeoDataFrame(pd.concat([canals, pipelines]), crs=f"EPSG:{epsg_jamaica}")

    # provide id
    nodes["node_id"] = nodes.apply(lambda node: f"{node.asset_type}_{node.OBJECTID}")
    edges["edge_id"] = edges.apply(lambda edge: f"{edge.asset_type}_{edge.OBJECTID}")

    # export as gpkg
    fname = os.path.join(
        processed_data_path, "networks", "water", "irrigation_assets_NIC.gpkg"
    )
    nodes.to_file(fname, layer="nodes", driver="GPKG")
    edges.to_file(fname, layer="edges", driver="GPKG")


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
