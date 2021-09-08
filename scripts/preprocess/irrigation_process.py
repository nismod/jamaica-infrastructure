"""Create the irrigation asset data for Jamaica
    Add the final asset data into a geopackage
"""
import sys
import os

import pandas as pd
import geopandas as gpd
import numpy as np
from preprocess_utils import *
from tqdm import tqdm
tqdm.pandas()

def main(config):
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    epsg_jamaica = 3448
    irrigation_data_path = os.path.join(incoming_data_path,'water','irrigation','raw')

    #### assumption list
    ## assume that reservoirs and dams (surface water sources) are not suscetible to flooding or hurricanes
    ## assume that the construction cost of potable water well is the same as irrigation well 
    ## assume that the flood damage function for canals is the same as for irrigation pipelines 

    # read in raw files, assign assign asset name for matching with cost/damage data, concat
    reservoirs = gpd.read_file(os.path.join(irrigation_data_path,'WATER SOURCES','RESERVOIRS.shp'))
    dams = gpd.read_file(os.path.join(irrigation_data_path,'WATER SOURCES','DAM.shp'))
    micro_dams = gpd.read_file(os.path.join(irrigation_data_path,'WATER SOURCES','MICRO_DAMS.shp'))
    wells = gpd.read_file(os.path.join(irrigation_data_path,'Well Sites','NIC_WELL_SITES_2020.shp'))
    canals = gpd.read_file(os.path.join(irrigation_data_path,'canal_network_NIC.shp'))[['OBJECTID','geometry']]
    pipelines = gpd.read_file(os.path.join(irrigation_data_path,'pipelines_network_NIC.shp'))[['OBJECTID','geometry']]

    wells['OBJECTID'] = np.linspace(0,wells.shape[0],wells.shape[0])
    wells = wells[['OBJECTID','geometry']]
    wells['Type'] = ['well']*wells.shape[0]
    wells['asset_type_cost_data'] = ['well']*wells.shape[0]
    wells['asset_type_flood_damage'] = ['well']*wells.shape[0]
    wells['asset_type_hurricane_damage'] = ['na']*wells.shape[0]
    
    canals['Type'] = ['canal']*canals.shape[0]
    canals['asset_type_cost_data'] = ['irrigation_canal']*canals.shape[0]
    canals['asset_type_flood_damage'] = ['irrigation_canal']*canals.shape[0]
    canals['asset_type_hurricane_damage'] = ['irrigation_canal']*canals.shape[0]
    
    pipelines['Type'] = ['pipeline']*pipelines.shape[0]
    pipelines['asset_type_cost_data'] = ['irrigation_pipeline']*pipelines.shape[0]
    pipelines['asset_type_flood_damage'] = ['irrigation_canal']*pipelines.shape[0]
    pipelines['asset_type_hurricane_damage'] = ['na']*pipelines.shape[0]

    irrigation_assets_NIC = pd.concat([wells,canals,pipelines])

    # # provide id
    irrigation_assets_NIC['Node_ID'] = irrigation_assets_NIC['OBJECTID'].astype(str) + irrigation_assets_NIC['Type'].astype(str)

    # # export as gpkg
    irrigation_assets_NIC = gpd.GeoDataFrame(irrigation_assets_NIC, 
                                                crs=f"EPSG:{epsg_jamaica}",
                                                geometry=irrigation_assets_NIC['geometry'])
    irrigation_assets_NIC.to_file(os.path.join(processed_data_path,
                                                    'networks',
                                                    'water',
                                                    'irrigation_assets_NIC.gpkg'),driver='GPKG')

if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
