"""Create the wastewater asset data for Jamaica
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
    wastewater_data_path = os.path.join(incoming_data_path,'water','wastewater')
    #### assumption list
    ## note sewers and storage tanks are assumed to not be susceptible to flooding/hurricane hazards and are therefore not provided with a cost or damage function
    ## wastewater treatment plants are assigned a cost based on capacity: limited capacity data was supplemented by data scraped from NWC website
    ## the cost of construction of wastewater treatment plants and pumping stations are assumed to be the same as for water treatment plants
    ## hurricane induced damage for water treatment plants is assumed to be the same as that for wastewater pumping stations and treatment plants
    ## assume cost of construction for sumps is equivalent to potable water intake and flood and hurricane damage curves are equivalent to potable pumping unit and water treatment plant respectively

    # read in raw files
    waste_water_facilities_NWC = gpd.read_file(os.path.join(wastewater_data_path,'raw','waste_water_facilities_NWC.shp'))
    sewers = gpd.read_file(os.path.join(wastewater_data_path,'raw','wGravityMain.shp'))

    # convert capacity data from gd to mgd
    waste_water_facilities_NWC['Capacity'] = pd.to_numeric(waste_water_facilities_NWC['Capacity'],errors='coerce')
    waste_water_facilities_NWC['capacity (mgd)'] = np.where(waste_water_facilities_NWC['Type']=='WW Treatment Plant', 
                                                    waste_water_facilities_NWC['Capacity']*0.000001, 
                                                    waste_water_facilities_NWC['Capacity'])

    # incorporate supplementary capacity data
    supplement_treatment_plant_capacity = pd.read_csv(os.path.join(wastewater_data_path,
                                                'additional_capacity_data',
                                                'supplement_treatment_plant_capacity.csv')) # merge on OBJECTID and TYPE
    waste_water_facilities_NWC = pd.merge(waste_water_facilities_NWC, supplement_treatment_plant_capacity, on='Name')

    # assign asset name for matching with cost/damage data
    asset_type_name_conversion = pd.DataFrame({'Type':['WW Relift Station', 'WW Pump Station', 'WW Treatment Plant', 'Sump'],
                                'asset_type_cost_data':['pumping unit','pumping unit','treatment plant','intake'],
                                'asset_type_flood_damage':['pumping unit','pumping unit','wastewater treatment plant','pumping unit'],
                                'asset_type_hurricane_damage':['water treatment plant','water treatment plant',
                                                                'water treatment plant','water treatment plant']})
    waste_water_facilities_NWC = pd.merge(waste_water_facilities_NWC, asset_type_name_conversion, on='Type')

    # provide id
    waste_water_facilities_NWC['Node_ID'] = waste_water_facilities_NWC['OBJECTID'].astype(str) + waste_water_facilities_NWC['Type'].astype(str)

    # export as gpkg
    waste_water_facilities_NWC = gpd.GeoDataFrame(waste_water_facilities_NWC, 
                                                    crs=f"EPSG:{epsg_jamaica}",
                                                    geometry=waste_water_facilities_NWC['geometry'])
    waste_water_facilities_NWC.to_file(os.path.join(processed_data_path,
                                                    'networks',
                                                    'water',
                                                    'waste_water_facilities_NWC.gpkg'),driver='GPKG')

if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)
