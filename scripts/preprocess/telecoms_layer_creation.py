"""Assign populations to electricity sinks for Jamaica
"""
import sys
import os

import pandas as pd
import geopandas as gpd
import numpy as np
import snkit
from preprocess_utils import *
from tqdm import tqdm
tqdm.pandas()

def match_areas_and_add_values(gdf_assign,gdf_values,gdf_assign_ids,gdf_values_ids,gdf_density_column):
    matches = gpd.sjoin(gdf_assign,
                        gdf_values, 
                        how="inner", predicate='intersects').reset_index()
    matches.rename(columns={"geometry":"assign_geometry"},inplace=True)
    matches = pd.merge(matches, gdf_values[gdf_values_ids+['geometry']],how="left",on=gdf_values_ids)
    matches["area_match"] = matches.progress_apply(lambda x:x["assign_geometry"].intersection(x["geometry"].buffer(0)).area,axis=1)
    matches[gdf_density_column] = matches["area_match"]*matches[gdf_density_column]

    matches = matches.groupby(gdf_assign_ids)[gdf_density_column].sum().reset_index()
    return pd.merge(gdf_assign,matches,how="left",on=gdf_assign_ids)



def main(config):
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    epsg_jamaica = 3448
    telecoms_nodes = gpd.read_file(os.path.join(incoming_data_path,
                                    "nsdmb", 
                                    "GWP_Jamaica_NSP_Master_Geodatabase_v01.gdb"), 
                        layer="cell_towers")
    telecoms_nodes = telecoms_nodes.to_crs(epsg=epsg_jamaica)
    telecoms_nodes["node_id"] = telecoms_nodes.index.values.tolist()
    telecoms_nodes["node_id"] = telecoms_nodes.progress_apply(lambda x: f"teln_{x.node_id}",axis=1)
    telecoms_nodes = snkit.network.drop_duplicate_geometries(telecoms_nodes, keep='first')
    telecoms_nodes = telecoms_nodes.reset_index()

    telecoms_nodes.to_file(os.path.join(processed_data_path,"networks","telecoms","telecoms.gpkg"),
                        layer="nodes",driver="GPKG")
    del telecoms_nodes
    telecoms_voronoi = gpd.read_file(os.path.join(incoming_data_path,"telecoms","telecoms_voronoi.gpkg"))
    population_areas = gpd.read_file(os.path.join(processed_data_path,
                            'population',
                            'population_projections.gpkg'),layer="mean")
    telecoms_voronoi = telecoms_voronoi.to_crs(epsg=epsg_jamaica)
    population_areas = population_areas.to_crs(epsg=epsg_jamaica)
    population_areas["2019"] = population_areas["2019"]/population_areas["AREA"]
    telecoms_voronoi = match_areas_and_add_values(telecoms_voronoi,population_areas,["node_id"],["ED_ID","ED"],"2019")
    telecoms_voronoi.to_file(os.path.join(processed_data_path,"networks","telecoms","telecoms.gpkg"),
                        layer="voronoi",driver="GPKG")


if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)