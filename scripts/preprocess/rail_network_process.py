"""Create the rail network costs for Jamaica
    The costs are derived from a World Bank PPI database
    We have extracted projects in the LAC region only
    We estimate the min-mean-max values from the database
    
    Write the costs is US$/km to the rail geopackage layer

"""
import sys
import os

import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt

from preprocess_utils import *
from tqdm import tqdm
tqdm.pandas()

def identify_rehabilitation_projects(x):
    if 'rehabilitate' in str(x['stype']).lower():
        return 1
    else:
        return 0


def main(config):
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    epsg_jamaica = 3448

    ppi_data = pd.read_stata(os.path.join(incoming_data_path,'world_bank_ppi','2019-Stata-Annual-sum.dta'))
    ppi_data = ppi_data[(ppi_data['Region'] == 'LAC') & (
                        ppi_data['income'] == 'Upper middle income') & (
                        ppi_data['ssector'] == 'Railroads') & (
                        ppi_data['status_n'] == 'Active')]
    ppi_data['rehab_tag'] = ppi_data.progress_apply(lambda x: identify_rehabilitation_projects(x),axis=1)
    ppi_data = ppi_data[ppi_data['rehab_tag'] == 1]
    ppi_data = ppi_data[(ppi_data['capacity'] == 'KM') & (ppi_data['pcapacity'] > 0) & (ppi_data['realphysicalassets'] > 0)]
    ppi_data = ppi_data.groupby('pcapacity')['realphysicalassets'].sum().reset_index()
    ppi_data['cost_per_km'] = ppi_data['realphysicalassets']/ppi_data['pcapacity']

    min_cost = 1e6*ppi_data.cost_per_km.quantile(0.10)
    average_cost = 1e6*ppi_data.cost_per_km.quantile(0.50)
    max_cost = 1e6*ppi_data.cost_per_km.quantile(0.90)

    edges = gpd.read_file(os.path.join(processed_data_path,
                                    'networks',
                                    'transport',
                                    'rail.gpkg'),
                            layer='edges')
    edges = gpd.GeoDataFrame(edges,geometry='geometry',crs={'init': f'epsg:{epsg_jamaica}'})
    edges['rail_length_m'] = edges.progress_apply(lambda x:x.geometry.length,axis=1)
    edges['cost_unit'] = 'USD/km'
    edges['min_damage_cost'] = min_cost
    edges['mean_damage_cost'] = average_cost
    edges['max_damage_cost'] = max_cost
    edges.to_file(os.path.join(processed_data_path,
                                    'networks',
                                    'transport',
                                    'rail.gpkg'),
                            layer='edges',driver="GPKG")

if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)