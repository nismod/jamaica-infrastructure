"""Map mining flows onto the rail network of Jamaica
"""

import sys
import os

import pandas as pd
import geopandas as gpd
import numpy as np
from shapely.geometry import LineString
from utils import *
from tqdm import tqdm

tqdm.pandas()
epsg_jamaica = 3448


def main(config):
    incoming_data_path = config["paths"]["incoming_data"]
    processed_data_path = config["paths"]["data"]

    airports = gpd.read_file(
        os.path.join(
            processed_data_path, "networks", "transport", "airport_polygon.gpkg"
        ),
        layer="areas",
    )
    airports = airports[airports["asset_type"] == "terminal"]
    airports = polygon_to_points(airports)
    airports["mode"] = "air"

    ports = gpd.read_file(
        os.path.join(processed_data_path, "networks", "transport", "port_polygon.gpkg"),
        layer="areas",
    )
    ports = polygon_to_points(ports)
    ports["mode"] = "port"

    rail_nodes = gpd.read_file(
        os.path.join(processed_data_path, "networks", "transport", "rail.gpkg"),
        layer="nodes",
    )
    # rail_nodes.loc[rail_nodes["node_id"]=="railn_124","asset_type"] = "station"
    # rail_nodes.loc[rail_nodes["node_id"]=="railn_44","status"] = "Functional"
    rail_nodes = rail_nodes[
        (rail_nodes["asset_type"] == "station") & (rail_nodes["status"] == "Functional")
    ]
    rail_nodes = rail_nodes.to_crs(epsg=epsg_jamaica)
    rail_nodes["mode"] = "rail"

    road_nodes = gpd.read_file(
        os.path.join(processed_data_path, "networks", "transport", "roads.gpkg"),
        layer="nodes",
    )
    road_nodes = road_nodes[road_nodes["component_id"] == 1]
    road_nodes = road_nodes.to_crs(epsg=epsg_jamaica)
    road_nodes["mode"] = "road"

    """Creating linkages
    """
    multi_modal = []
    multi_modal.append(
        map_nearest_locations_and_create_lines(
            airports.copy(), road_nodes.copy(), "node_id", "node_id", "air", "road"
        )
    )
    multi_modal.append(
        map_nearest_locations_and_create_lines(
            ports.copy(), road_nodes.copy(), "node_id", "node_id", "port", "road"
        )
    )
    multi_modal.append(
        map_nearest_locations_and_create_lines(
            ports.copy(), rail_nodes.copy(), "node_id", "node_id", "port", "rail"
        )
    )
    multi_modal.append(
        map_nearest_locations_and_create_lines(
            rail_nodes.copy(), road_nodes.copy(), "node_id", "node_id", "rail", "road"
        )
    )

    """Add road and rail
    """

    multi_modal = gpd.GeoDataFrame(
        pd.concat(multi_modal, axis=0, ignore_index=True),
        geometry="geometry",
        crs=f"EPSG:{epsg_jamaica}",
    )
    multi_modal["edge_id"] = multi_modal.index.values.tolist()
    multi_modal["edge_id"] = multi_modal.progress_apply(
        lambda x: f"multie_{x.edge_id}", axis=1
    )
    multi_modal["speed"] = 40.0
    multi_modal["time"] = 0.001 * multi_modal["length_m"] / multi_modal["speed"]
    print(multi_modal)

    rail_edges = gpd.read_file(
        os.path.join(processed_data_path, "networks", "transport", "rail.gpkg"),
        layer="edges",
    )
    rail_edges = rail_edges[rail_edges["status"] == "Functional"]
    rail_edges["from_mode"] = "rail"
    rail_edges["to_mode"] = "rail"
    rail_edges["time"] = 0.001 * rail_edges["length_m"] / rail_edges["speed"]

    road_edges = gpd.read_file(
        os.path.join(processed_data_path, "networks", "transport", "roads.gpkg"),
        layer="edges",
    )
    road_edges_max_id = max(
        [int(c.split("_")[1]) for c in road_edges["edge_id"].values.tolist()]
    )
    road_edges = road_edges[road_edges["component_id"] == 1]
    road_edges["from_mode"] = "road"
    road_edges["to_mode"] = "road"
    """Add road to fix lack of connectivity
    """

    road_edges["time"] = 0.001 * road_edges["length_m"] / road_edges["speed_kph"]
    road_edges = road_edges.rename(columns={"speed_kph": "speed"})
    print(rail_edges)
    print(road_edges)
    columns = [
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
    multi_modal = gpd.GeoDataFrame(
        pd.concat(
            [
                multi_modal[columns],
                rail_edges[columns],
                road_edges[columns],
            ],
            axis=0,
            ignore_index=True,
        ),
        geometry="geometry",
        crs=f"EPSG:{epsg_jamaica}",
    )
    print(multi_modal)
    multi_modal.to_file(
        os.path.join(
            processed_data_path, "networks", "transport", "multi_modal_network.gpkg"
        ),
        layer="edges",
        driver="GPKG",
    )
    columns = ["node_id", "mode", "geometry"]
    multi_modal = gpd.GeoDataFrame(
        pd.concat(
            [
                airports[columns],
                ports[columns],
                rail_nodes[columns],
                road_nodes[columns],
            ],
            axis=0,
            ignore_index=True,
        ),
        geometry="geometry",
        crs=f"EPSG:{epsg_jamaica}",
    )
    print(multi_modal)
    multi_modal.to_file(
        os.path.join(
            processed_data_path, "networks", "transport", "multi_modal_network.gpkg"
        ),
        layer="nodes",
        driver="GPKG",
    )


if __name__ == "__main__":
    CONFIG = load_config()
    main(CONFIG)
