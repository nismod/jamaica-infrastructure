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
    epsg_jamaica = 3448

    road_edges_dir = os.path.join(database_path,'National Land Agency Datasets','PIOJ','Street Centerline')
    roads = []
    for root, dirs, files in os.walk(road_edges_dir):
        for file in files:
            if file.endswith(".shp"):
                roads.append(gpd.read_file(os.path.join(road_edges_dir,file)))

    roads = gpd.GeoDataFrame(pd.concat(roads,ignore_index=True,sort='False',axis=0),
                            geometry='geometry',crs={'init': f'epsg:{epsg_jamaica}'})            
    
    roads = roads[~roads['CLASS'].isna()]
    # print (roads)
    roads = roads[~roads['geometry'].isna()]
    roads.to_file(os.path.join(os.path.join(database_path,'pioj_roads','detailed_roads.shp')))           
    # print (roads) 

    bridges = gpd.read_file(os.path.join(database_path,'NWA Bridges 2016','NWA Bridges 2016.shp'))
    bridges = bridges.to_crs(epsg=epsg_jamaica)

    create_network_from_nodes_and_edges(bridges,roads,
                                'road',
                                os.path.join(processed_data_path,'networks','transport','detailed_roads.gpkg'))

    

if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)