"""Assign populations to electricity sinks for Jamaica
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

    electricity_nodes = gpd.read_file(os.path.join(processed_data_path,
                            'networks','energy',
                            'electricity_network_v0.1.gpkg'),layer='nodes_processed')
    electricity_nodes = electricity_nodes[electricity_nodes['asset_type'] == 'sink']
    electricity_id_column = 'id'


    parish_boundaries = gpd.read_file(os.path.join(processed_data_path,
                            'boundaries',
                            'admin_boundaries.gpkg'),layer='admin1')
    
    population_boundaries = gpd.read_file(os.path.join(processed_data_path,
                            'boundaries',
                            'admin_boundaries.gpkg'),layer='admin3')
    # There is a population column = TOTAL_POP
    population_value_column = 'TOTAL_POP' 

    # For the population column = in the layer several entries of this layer are 0, which might not be right
    # There are also two other columns TOTAL_FMLE & TOTAL_MLE, which we might to add to get better estimates
    # population_boundaries['population'] = population_boundaries['TOTAL_FMLE'] + population_boundaries['TOTAL_MLE']
    # population_value_column = 'population'

    electricity_nodes = electricity_nodes.to_crs(epsg=epsg_jamaica)
    population_boundaries = population_boundaries.to_crs(epsg=epsg_jamaica)
    """We just take the whole electricity sinks and population 
    """
    electricity_nodes = assign_node_weights_by_population_proximity(electricity_nodes,
                        population_boundaries,
                        electricity_id_column,population_value_column,epsg=epsg_jamaica)
    print (electricity_nodes)

if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)