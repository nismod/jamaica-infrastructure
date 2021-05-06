"""Functions for preprocessing data
"""
import sys
import os
import json

import pandas as pd
import geopandas as gpd
import fiona
import numpy as np
import snkit

def load_config():
    """Read config.json
    """
    config_path = os.path.join(os.path.dirname(__file__), '..', '..', 'config.json')
    with open(config_path, 'r') as config_fh:
        config = json.load(config_fh)
    return config

def geopandas_read_file_type(file_path,file_layer,file_database=None):
    if file_database is not None:
        return gpd.read_file(os.path.join(file_path,file_database),
                        layer=file_layer)
    else:
        return gpd.read_file(os.path.join(file_path,file_layer))


def create_network_from_nodes_and_edges(nodes,edges,node_edge_prefix,out_fname):
    edges.columns = map(str.lower,  edges.columns)
    if 'id' in edges.columns.values.tolist():
        edges.rename(columns={'id':'e_id'},inplace=True)

    network = snkit.Network(nodes, edges)
    print ('* Done with network creation')

    network = snkit.network.split_multilinestrings(network)
    print ('* Done with splitting multilines')

    network = snkit.network.snap_nodes(network)
    print ('* Done with snapping nodes to edges')

    network = snkit.network.split_edges_at_nodes(network)
    print ('* Done with splitting edges at nodes')

    network = snkit.network.add_endpoints(network)   
    print ('* Done with adding endpoints')

    network = snkit.network.add_ids(network, 
                            edge_prefix=f"{node_edge_prefix}e", 
                            node_prefix=f"{node_edge_prefix}n")
    network = snkit.network.add_topology(network, id_col='id')
    print ('* Done with network topology')

    network.edges.rename(columns={'from_id':'from_node',
                                'to_id':'to_node',
                                'id':'edge_id'},
                                inplace=True)
    network.nodes.rename(columns={'id':'node_id'},inplace=True)
    
    network.edges.to_file(out_fname, layer='edges', driver='GPKG')
    network.nodes.to_file(out_fname, layer='nodes', driver='GPKG')

    return network
