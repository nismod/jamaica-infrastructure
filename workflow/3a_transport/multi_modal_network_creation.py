"""
Connect ports, airports, railways and roads together into a multi-modal transport network.
"""

import logging

import click
import pandas as pd
import geopandas as gpd
from tqdm import tqdm

from jamaica_infrastructure.transport.utils import (
    map_nearest_locations_and_create_lines,
    polygon_to_points,
)

tqdm.pandas()
epsg_jamaica = 3448


@click.command()
@click.version_option("1.0.0")
@click.option(
    "--road-network-path", "-r", required=True,
    type=click.Path(dir_okay=False, file_okay=True, exists=True),
    help="Path to road network geopackage.",
)
@click.option(
    "--rail-network-path", "-l", required=True,
    type=click.Path(dir_okay=False, file_okay=True, exists=True),
    help="Path to rail network geopackage.",
)
@click.option(
    "--port-path", "-p", required=True,
    type=click.Path(dir_okay=False, file_okay=True, exists=True),
    help="Path to ports geopackage.",
)
@click.option(
    "--airport-path", "-a", required=True,
    type=click.Path(dir_okay=False, file_okay=True, exists=True),
    help="Path to airports geopackage.",
)
@click.option(
    "--output-path", "-o", required=True,
    type=click.Path(dir_okay=False, file_okay=True, exists=False),
    help="Path to write multi-modal network geopackage.",
)
def main(
    road_network_path: str,
    rail_network_path: str,
    port_path: str,
    airport_path: str,
    output_path: str
):

    airports = gpd.read_file(airport_path, layer="areas")
    airports = airports[airports["asset_type"] == "terminal"]
    airports = polygon_to_points(airports)
    airports["mode"] = "air"

    ports = gpd.read_file(port_path, layer="areas")
    ports = polygon_to_points(ports)
    ports["mode"] = "port"

    rail_nodes = gpd.read_file(rail_network_path, layer="nodes")
    rail_nodes = rail_nodes[
        (rail_nodes["asset_type"] == "station") & (rail_nodes["status"] == "Functional")
    ]
    rail_nodes = rail_nodes.to_crs(epsg=epsg_jamaica)
    rail_nodes["mode"] = "rail"

    road_nodes = gpd.read_file(road_network_path, layer="nodes")
    road_nodes = road_nodes[road_nodes["component_id"] == 1]
    road_nodes = road_nodes.to_crs(epsg=epsg_jamaica)
    road_nodes["mode"] = "road"

    # create linkages
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

    # road and rail
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

    rail_edges = gpd.read_file(rail_network_path, layer="edges")
    rail_edges = rail_edges[rail_edges["status"] == "Functional"]
    rail_edges["from_mode"] = "rail"
    rail_edges["to_mode"] = "rail"
    rail_edges["time"] = 0.001 * rail_edges["length_m"] / rail_edges["speed"]

    road_edges = gpd.read_file(road_network_path, layer="edges")
    road_edges = road_edges[road_edges["component_id"] == 1]
    road_edges["from_mode"] = "road"
    road_edges["to_mode"] = "road"

    road_edges["time"] = 0.001 * road_edges["length_m"] / road_edges["speed_kph"]
    road_edges = road_edges.rename(columns={"speed_kph": "speed"})

    logging.info(f"Writing multi-modal network to disk: {output_path}")
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
    multi_modal.to_file(output_path, layer="edges", driver="GPKG")

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
    multi_modal.to_file(output_path, layer="nodes", driver="GPKG")


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    main()
