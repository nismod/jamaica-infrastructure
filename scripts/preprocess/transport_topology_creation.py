"""Take road and rail data from the NSDMB database and create topological networks
    Add these networks into a geopackage
"""
import sys
import os

import geopandas as gpd
from preprocess_utils import *


def main(config):
    database_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    database_name = 'GWP_Jamaica_NSP_Master_Geodatabase_v01.gdb'

    sector_details = [
                        {
                            'sector_name':'rail',
                            'node_folder':'nsdmb',
                            'node_layer':'railway_stations',
                            'node_file_type':database_name,
                            'edge_folder':'nsdmb',
                            'edge_layer':'railways',
                            'edge_file_type':database_name,

                        },
                        {
                            'sector_name':'road',
                            'node_folder':'NWA Bridges 2016',
                            'node_layer':'NWA Bridges 2016.shp',
                            'node_file_type':None,
                            'edge_folder':'nsdmb',
                            'edge_layer':'roads_main_NWA',
                            'edge_file_type':database_name,

                        }

                    ]

    for sector in sector_details:
        if sector['node_layer']:
            nodes = geopandas_read_file_type(os.path.join(database_path,
                                                sector['node_folder']
                                                ),
                                            sector['node_layer'],
                                            sector['node_file_type'])
        else:
            nodes = None
        if sector['edge_layer']:
            edges = geopandas_read_file_type(os.path.join(database_path,
                                                sector['edge_folder']
                                                ),
                                            sector['edge_layer'],
                                            sector['edge_file_type'])
        else:
            print ('{} error: there should be at an edge layer'.format(sector['sector_name']))
            break

        create_network_from_nodes_and_edges(nodes,edges,
                                    sector['sector_name'],
                                    os.path.join(processed_data_path,'networks',f"{sector['sector_name']}.gpkg"))

    

if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)