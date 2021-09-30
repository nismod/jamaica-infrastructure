"""Create the pipelines asset data for Jamaica
    Add the final asset data into a geopackage
"""
import os

import pandas as pd
import geopandas as gpd
import numpy as np
from preprocess_utils import *
from tqdm import tqdm
tqdm.pandas()

def main(config):
    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]
    epsg_jamaica = 3448
    potable_data_path = os.path.join(
        incoming_data_path, "water", "potable", "raw"
    )
    sewer_data_path = os.path.join(
        incoming_data_path, "water", "wastewater", "raw"
    )
    #### assumption list
    # potable pipes and sewers are same cost

    # # read in raw files
    pot_pipelines = gpd.read_file(os.path.join(potable_data_path,'pipelines_network_NWC.shp'))[
                                            ['OBJECTID','Material','Diam_Nom','Shape_Leng','geometry']]
    sewer_g = gpd.read_file(os.path.join(sewer_data_path,'wGravityMain.shp'))[
                        ['OBJECTID','Material','Nom_Diam','Shape_Leng','geometry']]
    sewer_p = gpd.read_file(os.path.join(sewer_data_path,'wPressurizedMain.shp'))[
                        ['OBJECTID','Material','Diam_Nom','Shape_Leng','geometry']]
    cost_data = pd.read_csv(os.path.join(incoming_data_path,"water","cost","water_asset_costs.csv"))

    pot_pipelines['Type'] = ['potable']*pot_pipelines.shape[0]
    pot_pipelines.columns = ['OBJECTID','Material','Diameter','Length','geometry','Type']
    sewer_g['Type'] = ['sewer_gravity']*sewer_g.shape[0]
    sewer_g.columns = ['OBJECTID','Material','Diameter','Length','geometry','Type']
    sewer_p['Type'] = ['sewer_pressure']*sewer_p.shape[0]
    sewer_p.columns = ['OBJECTID','Material','Diameter','Length','geometry','Type']
    pipelines = pd.concat([pot_pipelines,sewer_g,sewer_p],axis=0,ignore_index=True)
    pipelines['Diameter'] = np.where(pipelines['Diameter']<0, 2, pipelines['Diameter'])
    
    # # assign asset name for matching with cost/damage data

    pipelines['asset_type_cost_data'] = ['pipeline']*pipelines.shape[0]

    pipelines = pd.merge(pipelines, cost_data, left_on = 'asset_type_cost_data', right_on = 'asset', how='left')
    # print (pipelines[['cost ($J) - lower bound','cost ($J) - upper bound','Diameter','Length']])
    
    pipelines['min_damage_cost'] = pipelines['cost ($J) - lower bound'].str.replace(",","").astype(float)*pipelines['Diameter']
    pipelines['max_damage_cost'] = pipelines['cost ($J) - upper bound'].str.replace(",","").astype(float)*pipelines['Diameter']
    pipelines.rename(columns={"Type": "asset_type"},inplace=True)
    pipelines['cost_unit'] = '$J/m'

    # provide id
    pipelines["edge_id"] = pipelines \
        .apply(lambda edge: f"{edge.asset_type}_{edge.OBJECTID}", axis=1)
    # # export as gpkg
    pipelines = gpd.GeoDataFrame(pipelines, crs=f"EPSG:{epsg_jamaica}",geometry=pipelines['geometry'])
    pipelines = pipelines[pipelines['Length'] > 0]
    # export as gpkg
    fname = os.path.join(
        processed_data_path, "networks", "water", "pipelines_NWC.gpkg"
    )
    pipelines.to_file(fname,layer = 'edges',driver='GPKG')

if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
