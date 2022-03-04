#!/usr/bin/env python
# coding: utf-8
"""Clean duplicate network elements from preprocessed data

Assume this script is run from processed_data/networks as current working 
directory
"""
import os
import warnings

import geopandas
import pandas


warnings.filterwarnings(action="ignore", module="geopandas")

pandas.read_csv("network_layers_hazard_intersections_details.csv")
files = pandas.read_csv(
    "network_layers_hazard_intersections_details.csv",
    usecols=["path", "asset_layer", "asset_id_column"],
)

# Report any duplicate row IDs
for _, layer, id_col, path in files.itertuples():
    if "buildings" in path:
        continue
    df = geopandas.read_file(os.path.join("..", path), layer=layer)
    count = len(df[id_col])
    unique = len(df[id_col].unique())
    if count != unique:
        print(layer, path, id_col, count, unique)


def find_dups(df, id_col):
    return df[df.duplicated(subset=[id_col], keep=False)]

# Drop irrigation edges
# NOTE drops columns with relation to nodes - should be in separate one- or 
# many-to-many mapping
water_irrigation_edges = (
    geopandas.read_file("water/irrigation_assets_NIC.gpkg", layer="edges")
    .drop(columns=["node_id", "sqm_per_node", "GDP/day_per_node", "OBJECTID"])
    .drop_duplicates(subset=["edge_id"])
    .reset_index(drop=True)
)
water_irrigation_edges.edge_id = water_irrigation_edges.edge_id.str.replace(
    ".0", "", regex=False
)
water_irrigation_edges.to_file("water/irrigation_assets_NIC.gpkg", layer="edges")

# Drop potable nodes
# NOTE population/GDP assignment columns seem to vary across duplicates - store mapping
# separately.
potable_nodes = (
    geopandas.read_file("water/potable_facilities_NWC.gpkg", layer="nodes")
    .drop(columns=["X", "Y", "lat", "lon", "OBJECTID", "LOCATION_y"])
    .rename(columns={"LOCATION_x": "LOCATION"})
    .drop_duplicates(subset=["node_id"])
    .reset_index(drop=True)
)
potable_nodes.node_id = (
    potable_nodes.node_id.str.lower()
    .str.replace(".0", "", regex=False)
    .str.replace(" ", "_", regex=False)
)
potable_nodes.to_file("water/potable_facilities_NWC.gpkg", layer="nodes")

# Drop electricity node
# NOTE single node was assigned to two parishes, otherwise identical
elec_nodes = geopandas.read_file(
    "energy/electricity_network_v3.0.gpkg", layer="nodes"
).drop_duplicates(subset=["id"])
elec_nodes.to_file("energy/electricity_network_v3.0.gpkg", layer="nodes")
