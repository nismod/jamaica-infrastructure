"""Create the irrigation asset data for Jamaica
    Add the final asset data into a geopackage
"""
import os

import pandas as pd
import geopandas as gpd
import numpy as np
from preprocess_utils import *

def estimate_costs_and_units(x):
    if x.asset_type == 'pipeline':
        min_damage_cost = float(str(x['cost ($J) - lower bound']).replace(",",""))/x.geometry.length
        max_damage_cost = float(str(x['cost ($J) - upper bound']).replace(",",""))/x.geometry.length
        cost_unit = '$J/m'
    else:
        min_damage_cost = float(str(x['cost ($J) - lower bound']).replace(",",""))
        max_damage_cost = float(str(x['cost ($J) - upper bound']).replace(",",""))
        cost_unit = '$J'

    return cost_unit,min_damage_cost,max_damage_cost



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

    cost_data = pd.read_csv(
        os.path.join(
            incoming_data_path,"water","cost","water_asset_costs.csv"
        )
    )

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

    nodes = gpd.GeoDataFrame(wells, crs=f"EPSG:{epsg_jamaica}",geometry='geometry')
    nodes = pd.merge(
        nodes, cost_data, 
        left_on = 'asset_type_cost_data',
        right_on = 'asset', how='left'
    )
    nodes['min_damage_cost'] = nodes['cost ($J) - lower bound'].str.replace(",","").astype(float)
    nodes['max_damage_cost'] = nodes['cost ($J) - upper bound'].str.replace(",","").astype(float)
    nodes['cost_unit'] = '$J'

    edges = pd.concat([canals, pipelines])
    edges = gpd.GeoDataFrame(edges, crs=f"EPSG:{epsg_jamaica}",geometry='geometry')
    edges = pd.merge(
        edges, cost_data, 
        left_on = 'asset_type_cost_data',
        right_on = 'asset', how='left'
    )
    edges['cost_and_units'] = edges.progress_apply(lambda x: estimate_costs_and_units(x),axis=1)
    edges[['cost_unit','min_damage_cost','max_damage_cost']] = edges['cost_and_units'].apply(pd.Series)
    edges.drop('cost_and_units',axis=1,inplace=True)

    # provide id
    nodes["node_id"] = nodes.apply(lambda node: f"{node.asset_type}_{node.OBJECTID}", axis=1)
    edges["edge_id"] = edges.apply(lambda edge: f"{edge.asset_type}_{edge.OBJECTID}", axis=1)

    # export as gpkg
    fname = os.path.join(
        processed_data_path, "networks", "water", "irrigation_assets_NIC.gpkg"
    )
    nodes.to_file(fname, layer="nodes", driver="GPKG")
    edges.to_file(fname, layer="edges", driver="GPKG")


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
