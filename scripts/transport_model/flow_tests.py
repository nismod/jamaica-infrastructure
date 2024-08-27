"""Map mining flows onto the rail network of Jamaica
"""

import sys
import os

import pandas as pd
import geopandas as gpd
import numpy as np
import igraph as ig
from utils import *
from tqdm import tqdm

tqdm.pandas()
epsg_jamaica = 3448


def route_areas_to_nearest_ports(
    areas,
    areas_id,
    areas_gdp,
    nodes,
    edges,
    ports,
    connection_type="areas",
    flow_type="exports",
):
    nearest_roads = map_nearest_locations_and_create_lines(
        areas.copy(),
        nodes[nodes["mode"] == "road"].copy(),
        areas_id,
        "node_id",
        connection_type,
        "road",
    )
    nearest_roads["edge_id"] = nearest_roads.progress_apply(
        lambda x: f"{connection_type}roade_{x.name}", axis=1
    )
    nearest_roads["speed"] = 10.0
    nearest_roads["time"] = 0.001 * nearest_roads["length_m"] / nearest_roads["speed"]

    roads_columns = [
        "from_node",
        "to_node",
        "edge_id",
        "from_mode",
        "to_mode",
        "length_m",
        "speed",
        "time",
        "geometry",
    ]
    network = pd.concat(
        [edges[roads_columns], nearest_roads[roads_columns]], axis=0, ignore_index=True
    )[roads_columns]
    G = ig.Graph.TupleList(
        network.itertuples(index=False), edge_attrs=list(network.columns)[2:]
    )

    all_ports = ports["node_id"].values.tolist()

    # od_pairs = [list(zip([b]*len(all_ports),all_ports)) for b in areas[areas_id].values.tolist()]
    od_pairs = [
        list(zip(all_ports, [b] * len(all_ports)))
        for b in areas[areas_id].values.tolist()
    ]
    od_pairs = [item for sublist in od_pairs for item in sublist]
    od_pairs = pd.DataFrame(od_pairs, columns=["origin_id", "destination_id"])
    od_pairs = pd.merge(
        od_pairs,
        areas[[areas_id, areas_gdp]],
        how="left",
        left_on=["destination_id"],
        right_on=[areas_id],
    )
    flow_paths = network_od_paths_assembly(
        od_pairs[["origin_id", "destination_id", areas_gdp]], G, "time", areas_gdp
    )
    flow_paths = flow_paths.sort_values(by="gcost")
    flow_paths = flow_paths.drop_duplicates(subset=["destination_id"], keep="first")
    flow_paths = flow_paths[
        ["origin_id", "destination_id", "edge_path", "gcost", areas_gdp]
    ]
    if flow_type == "exports":
        flow_paths.columns = [
            "destination_id",
            "origin_id",
            "edge_path",
            "gcost",
            areas_gdp,
        ]

    return flow_paths, nearest_roads[roads_columns]


def main(config):
    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]
    results_path = config["paths"]["output"]

    """Get the ports and multimodal network
    """
    ports = gpd.read_file(
        os.path.join(processed_data_path, "networks", "transport", "port_polygon.gpkg"),
        layer="areas",
    )
    ports["export_tonnes"] = ports["export_tonnes"].fillna(0)
    ports = ports[(ports["commodity"] != "alumina") & (ports["export_tonnes"] > 0)]
    ports = ports.drop_duplicates(subset=["name"], keep="first")
    ports = ports[["node_id", "name", "export_tonnes", "geometry"]]
    ports["export_wt"] = ports["export_tonnes"] / ports["export_tonnes"].sum()
    ports["geometry"] = ports.progress_apply(lambda x: x.geometry.centroid, axis=1)
    ports = ports.to_crs(epsg=epsg_jamaica)
    print(ports)

    nodes = gpd.read_file(
        os.path.join(
            processed_data_path, "networks", "transport", "multi_modal_network.gpkg"
        ),
        layer="nodes",
    )
    nodes = nodes.to_crs(epsg=epsg_jamaica)
    edges = gpd.read_file(
        os.path.join(
            processed_data_path, "networks", "transport", "multi_modal_network.gpkg"
        ),
        layer="edges",
    )
    edges = edges[(edges["from_mode"] != "rail") & (edges["to_mode"] != "rail")]
    print(edges)

    source = "roadsn_84994"
    target = "roadsn_83456"
    target = "roadsn_33205"
    # target = "roadsn_84979"
    # target = "roadsn_83546"
    roads_columns = [
        "from_node",
        "to_node",
        "edge_id",
        "from_mode",
        "to_mode",
        "length_m",
        "speed",
        "time",
        "geometry",
    ]
    network = edges[roads_columns]
    G = ig.Graph.TupleList(
        network.itertuples(index=False), edge_attrs=list(network.columns)[2:]
    )
    road_path, road_cost = network_od_path_estimations(G, source, target, "time")
    print(road_path)


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
