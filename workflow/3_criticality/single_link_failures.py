"""Transport failure analysis with rerouting."""

import logging

import click
import pandas as pd
import geopandas as gpd
from tqdm import tqdm

from jamaica_infrastructure.transport.econ import economic_losses_from_network_damage
from jamaica_infrastructure.transport.flow import (
    igraph_scenario_edge_failures_premade_network,
    read_flow_data,
)


tqdm.pandas()


@click.command()
@click.version_option("1.0.0")
@click.option(
    "--edge-chunk-map-csv",
    "-c",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True),
    help="Path to the file mapping a chunk ID to start and stop edge indicies.",
)
@click.option("--chunk-id", "-i", required=True, type=int, help="Path to the file mapping a chunk ID to start and stop edge indicies.")
@click.option(
    "--edges-file",
    "-e",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Jamaican Multimodal Network GeoPackage file",
)
@click.option(
    "--flow-data-dir",
    "-f",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True),
    help="Path to the directory containing nominal flow data (ODs, etc.).",
)
@click.option(
    "--output-path",
    "-o",
    required=True,
    type=click.Path(exists=False, file_okay=True, dir_okay=False, readable=True),
    help="Path to write failure analysis results to.",
)
@click.option(
    "--hourly-wage",
    "-w",
    required=True,
    type=float,
    help="Hourly wage in JMD for labour cost calculations.",
)
@click.option(
    "--trade-effect",
    "-t",
    required=True,
    type=float,
    help="Fraction of trade value affected by rerouting.",
)
def main(*, edge_chunk_map_csv, chunk_id, edges_file, flow_data_dir, output_path, hourly_wage, trade_effect):
    """
    Calculate economic costs of removing road edges from multi-modal transport network.
    """

    chunk_map = pd.read_csv(edge_chunk_map_csv).set_index("id")
    min_edge_number, max_edge_number = sorted(chunk_map.loc[chunk_id])

    logging.info(f"Running transport failure analysis for edge positions {min_edge_number} -> {max_edge_number}")

    logging.info("Read nominal flow data")
    network_dictionary, all_flows, trade_sectors = read_flow_data(flow_data_dir)

    logging.info("Reading network data")
    edges = gpd.read_file(edges_file, layer="edges")

    rail_edges = edges[(edges["from_mode"] == "rail") & (edges["to_mode"] == "rail")]["edge_id"].values.tolist()
    road_edges = edges[(edges["from_mode"] == "road") & (edges["to_mode"] == "road")]["edge_id"].values.tolist()
    edge_fail = rail_edges + road_edges

    logging.info("Failing edges and reallocating flows")

    if max_edge_number > len(edge_fail):
        max_edge_number = len(edge_fail)

    edge_fail_results = []
    for edge_number in range(min_edge_number, max_edge_number):
        edge = edge_fail[edge_number]
        logging.info(f"Failing {edge}")
        for networks in network_dictionary.values():
            edge_fail_results += igraph_scenario_edge_failures_premade_network(
                # we will remove edges from the graph, only operate on a copy
                networks["network"].copy(),
                [edge],
                networks["flows"],
                networks["edge_indexes"],
                "edge_path",
                "time",
            )

    logging.info("Done failing edges")

    logging.info("Calculating resulting economic losses")
    edge_fail_results = pd.DataFrame(edge_fail_results)
    losses = economic_losses_from_network_damage(
        edge_fail_results,
        "edge_id",
        "no_access",
        all_flows,
        trade_sectors,
        trade_effect,
        hourly_wage
    )
    logging.info(f"Losses:\n{losses}")

    logging.info("Writing results to disk")
    losses.to_csv(output_path, index=False)


if __name__ == "__main__":

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)
    main()
