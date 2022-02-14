#!/usr/bin/env python
# coding: utf-8
"""Add connected component ID to network

Use this to sense-check whether a network is fully connected
"""
import sys
from glob import glob

import geopandas
import networkx


def main(network_fname, out_fname, node_id_col):
    edges = geopandas.read_file(network_fname, layer="edges")
    nodes = geopandas.read_file(network_fname, layer="nodes")

    G = networkx.Graph()
    G.add_nodes_from(
        (getattr(n, node_id_col), {"geometry": n.geometry}) for n in nodes.itertuples()
    )
    G.add_edges_from(
        (e.from_node, e.to_node, {"edge_id": e.edge_id, "geometry": e.geometry})
        for e in edges.itertuples()
    )
    components = networkx.connected_components(G)
    for num, c in enumerate(components):
        print(f"Component {num} has {len(c)} nodes")
        edges.loc[(edges.from_node.isin(c) | edges.to_node.isin(c)), "component"] = num
        nodes.loc[nodes[node_id_col].isin(c), "component"] = num

    edges.to_file(out_fname, layer="edges", driver="GPKG")
    nodes.to_file(out_fname, layer="nodes", driver="GPKG")

if __name__ == "__main__":
    network_fname, out_fname, node_id_col = sys.argv[1:]
    main(network_fname, out_fname, node_id_col)
