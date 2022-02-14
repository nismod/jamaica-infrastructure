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
import assign_network_components
from collections import defaultdict
from preprocess_utils import *
from tqdm import tqdm
tqdm.pandas()

def identify_rehabilitation_projects(x):
    if 'rehabilitate' in str(x['stype']).lower():
        return 1
    else:
        return 0

def get_node_status(x,edges):
    edge_status = list(set(edges[(edges.from_node == x.node_id) | (edges.to_node == x.node_id)]['status'].values.tolist()))
    if "Functional" in edge_status:
        return "Functional"
    else:
        return "Non-Functional"


def main(config):
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    epsg_jamaica = 3448

    database_name = "GWP_Jamaica_NSP_Master_Geodatabase_v01.gdb"

    edges = gpd.read_file(os.path.join(incoming_data_path,'nsdmb','GWP_Jamaica_NSP_Master_Geodatabase_v01.gdb'),
                        layer='railways')
    print (edges)
    nodes = gpd.read_file(os.path.join(incoming_data_path,'nsdmb','GWP_Jamaica_NSP_Master_Geodatabase_v01.gdb'),
                        layer='railway_stations')
    print (nodes)
    network = create_network_from_nodes_and_edges(
            nodes,
            edges,
            "rail"
        )
    print (network.edges)
    print (network.nodes)
    edges = network.edges
    nodes = network.nodes
    # edges = gpd.read_file(os.path.join(file_path, file_database), layer=file_layer)
    # average_cost = 0.001*945827.0/1.609 # Cost estimate obtained from JRC project report is 945,827 US$/mile
    min_cost = 0.001*1.0e6/1.609/0.0067 # Cost estimate obtained from Ministry is 1,000,000 US$/mile
    max_cost = 0.001*1.2*1.0e6/1.609/0.0067 # Cost estimate obtained from Ministry is 1,200,000 US$/mile 
    average_cost = 0.5*(min_cost + max_cost)
    # edges = gpd.read_file(os.path.join(processed_data_path,
    #                                 'networks',
    #                                 'transport',
    #                                 'rail.gpkg'),
    #                         layer='edges')
    edges = gpd.GeoDataFrame(edges,geometry='geometry',crs=f"EPSG:{epsg_jamaica}")
    edges.to_file(os.path.join(processed_data_path,
                                    'networks',
                                    'transport',
                                    'rail.gpkg'),
                            layer='edges',driver="GPKG")

    nodes = gpd.GeoDataFrame(nodes,geometry='geometry',crs=f"EPSG:{epsg_jamaica}")
    nodes['asset_type'] = nodes.progress_apply(lambda x: 'station' if str(x.Station) not in ("None","nan") else 'junction',axis=1)
    nodes['status'] = nodes.progress_apply(lambda x: get_node_status(x,edges),axis=1)
    nodes['cost_unit'] = 'J$'
    station_cost = 500000.0/0.0067 # Cost estimate obtained from JRC project report
    nodes['min_damage_cost'] = 0.8*station_cost
    nodes['mean_damage_cost'] = station_cost
    nodes['max_damage_cost'] = 1.2*station_cost
    nodes.to_file(os.path.join(processed_data_path,
                                    'networks',
                                    'transport',
                                    'rail.gpkg'),
                            layer='nodes',driver="GPKG")
    component_file = os.path.join(processed_data_path,
                                    'networks',
                                    'transport',
                                    'rail.gpkg')
    edges_components = assign_network_components.main(component_file,component_file,"node_id")

    nodes = gpd.read_file(os.path.join(processed_data_path,
                                    'networks',
                                    'transport',
                                    'rail.gpkg'),
                            layer='nodes')
    print (nodes)

    edges = gpd.read_file(os.path.join(processed_data_path,
                                    'networks',
                                    'transport',
                                    'rail.gpkg'),
                            layer='edges')
    max_id = max([int(v.split("_")[1]) for v in edges["edge_id"].values.tolist()])
    components = list(set(nodes["component"].values.tolist()))
    print (components)
    joined_components = []
    for i in range(len(components)-1): 
        node_component = nodes[nodes["component"] == components[i]][["node_id","geometry"]]
        remaining_components = nodes[nodes["component"].isin(components[i+1:])][["node_id","geometry"]]
        component_distances = map_nearest_locations_and_create_lines(node_component,remaining_components,"node_id","node_id","rail","rail")
        joined_components.append(component_distances[component_distances["length_m"] <= 150])

    joined_components = pd.concat(joined_components,axis=0,ignore_index=True)
    joined_components["edge_id"] = list(max_id + 1 + joined_components.index.values)
    joined_components["edge_id"] = joined_components.progress_apply(lambda x: f"raile_{x.edge_id}",axis=1)
    joined_components.drop(["from_mode","to_mode"],axis=1,inplace=True)
    print (joined_components)

    print (edges)
    joined_components["status"] = "Functional"
    for row in joined_components.itertuples():
        row_dict = defaultdict()
        row_dict["edge_id"] = row.edge_id
        from_values = nodes[nodes["node_id"] == row.from_node]["status"].values[0]
        to_values = nodes[nodes["node_id"] == row.to_node]["status"].values[0]
        values = list(set([from_values,to_values]))
        joined_components.loc[joined_components["edge_id"] == row.edge_id,"status"] = values[0]

    print (joined_components)

    edges = gpd.GeoDataFrame(pd.concat([edges,joined_components],axis=0,ignore_index=True),geometry="geometry",crs=f"EPSG:{epsg_jamaica}")
    edges['asset_type'] = edges.progress_apply(lambda x: 'rail' if x.status == 'Functional' else x.status, axis=1)
    edges['length_m'] = edges.progress_apply(lambda x:x.geometry.length,axis=1)
    edges['cost_unit'] = 'J$/m'
    edges['min_damage_cost'] = min_cost
    edges['mean_damage_cost'] = average_cost
    edges['max_damage_cost'] = max_cost
    edges['speed'] = 120 # Assumed speed for railways

    edges.to_file(os.path.join(processed_data_path,
                                    'networks',
                                    'transport',
                                    'rail.gpkg'),
                            layer='edges',driver="GPKG")

    nodes.loc[nodes["node_id"] == "railn_156","asset_type"] = "station"

    nodes.to_file(os.path.join(processed_data_path,
                                    'networks',
                                    'transport',
                                    'rail.gpkg'),
                            layer='nodes',driver="GPKG")
    component_file = os.path.join(processed_data_path,
                                    'networks',
                                    'transport',
                                    'rail.gpkg')
    edges_components = assign_network_components.main(component_file,component_file,"node_id")

if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)