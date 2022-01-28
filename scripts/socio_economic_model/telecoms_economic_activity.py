"""Assign agriculture GDP to land use layers in Jamaica
"""
import sys
import os
import subprocess

import rasterio
import rioxarray
import pandas as pd
import geopandas as gpd
# gpd._compat.USE_PYGEOS = True
# gpd.options.use_pygeos = True
from shapely.geometry import Point
import shapely
import numpy as np
from utils import *
from tqdm import tqdm
tqdm.pandas()

def match_buildings_to_footprints(buildings,gdf,building_id,gdf_id):
    matches = gpd.sjoin(
                        buildings,
                        gdf,
                        how="inner", predicate='intersects').reset_index()
    matches.rename(columns={"geometry":"building_geometry"},inplace=True)
    matches = pd.merge(matches, gdf[[gdf_id,'geometry']],how="left",on=gdf_id)
    matches["area_match"] = matches.progress_apply(lambda x:x["building_geometry"].intersection(x["geometry"].buffer(0)).area,axis=1)
    matches = matches.sort_values(by=["area_match"],ascending=False)
    matches = matches.drop_duplicates(subset=[building_id], keep="first")
    matches.drop(["area_match","building_geometry","geometry"],axis=1,inplace=True)

    return matches

def main(config):
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']

    epsg_jamaica = 3448
    post_office_share = 0.157
    financial_year = 2019
    population_year = 2019
    sector_codes = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O']
    gdp_columns = [f"{scode}_GDP" for scode in sector_codes]
    telecoms = gpd.read_file(os.path.join(processed_data_path,
                        "networks","telecoms",
                        "telecoms.gpkg"),
                        layer="voronoi")
    telecoms = telecoms.to_crs(epsg=epsg_jamaica)
    telecoms_id = "node_id"
    
    buildings = gpd.read_file(os.path.join(
                                    processed_data_path,
                                    "buildings",
                                    'buildings_assigned_economic_activity.gpkg'),
                                  layer="areas")
    buildings = buildings.to_crs(epsg=epsg_jamaica)
    building_id = "osm_id"

    economic_output = pd.read_excel(os.path.join(processed_data_path,
                                        "macroeconomic_data",
                                        "detailed_sector_GVA_GDP_current_prices.xlsx"),
                        sheet_name="2019")
    economic_output.columns = [str(c).strip() for c in economic_output.columns.values.tolist()]
    economic_output["subsector_code"] = economic_output["subsector_code"].apply(str)
    totat_gva = economic_output[economic_output["sector_code"] == "GVA"][f"{financial_year}"].sum()
    total_tax =  economic_output[economic_output["sector_code"] == "TAX"][f"{financial_year}"].sum()
    tax_rate = 1.0*total_tax/totat_gva

    output = (1.0 - post_office_share)*(1+tax_rate)*(1.0e6/365.0)*economic_output[(
                                        economic_output["sector_code"] == "I"
                                        ) & (
                                        economic_output["subsector_code"] == "640"
                                        )][f"{financial_year}"].sum()
    
    telecoms["T_GDP"] = output*telecoms[f"{population_year}"]/telecoms[f"{population_year}"].sum()
    telecoms["GDP_unit"] = "JD/day"
    telecoms.to_file(os.path.join(processed_data_path,
                        "networks","telecoms",
                        "telecoms.gpkg"),
                        layer="voronoi")

    telecoms_nodes = gpd.read_file(os.path.join(processed_data_path,
                        "networks","telecoms",
                        "telecoms.gpkg"),
                        layer="nodes")
    telecoms_nodes = pd.merge(telecoms_nodes,telecoms[[telecoms_id,"T_GDP"]],how="left",on=[telecoms_id])
    telecoms_nodes["T_GDP"] = telecoms_nodes["T_GDP"].fillna(0)
    telecoms_nodes["GDP_unit"] = "JD/day"
    telecoms_nodes.to_file(os.path.join(processed_data_path,
                        "networks","telecoms",
                        "telecoms.gpkg"),
                        layer="nodes")

    telecoms_buildings = match_buildings_to_footprints(buildings,telecoms[[telecoms_id,"geometry"]],building_id,telecoms_id)
    telecoms_buildings = telecoms_buildings[[telecoms_id,building_id] + gdp_columns]
    telecoms_buildings.to_csv(os.path.join(processed_data_path,
                                    "networks_economic_activity",
                                    "telecoms_buildings_economic_activity_mapping.csv"),
                            index=False)

    telecoms_buildings = telecoms_buildings.groupby([telecoms_id])[gdp_columns].sum().reset_index()
    telecoms_buildings = pd.merge(telecoms[[telecoms_id,"T_GDP"]],telecoms_buildings,how="left",on=[telecoms_id]).fillna(0)
    telecoms_buildings["total_GDP"] = telecoms_buildings[["T_GDP"] + gdp_columns].sum(axis=1)
    telecoms_buildings["GDP_unit"] = "JD/day"
    telecoms_buildings.to_csv(os.path.join(processed_data_path,
                                    "networks_economic_activity",
                                    "telecoms_dependent_economic_activity.csv"),
                            index=False)    


if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)