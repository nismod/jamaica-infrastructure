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
    "--edges-file",
    "-e",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Jamaican Multimodal Network GeoPackage file",
)
@click.option(
    "--rail-nodes-file",
    "-r",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Railway nodes GeoPackage file",
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
def main(*, flow_data_dir, edges_file, rail_nodes_file, output_path, hourly_wage, trade_effect):
    """
    Calculate economic costs of removing railway stations and adjacent edges
    from multi-modal transport network.
    """

    logging.info("Read nominal flow data")
    network_dictionary, all_flows, trade_sectors = read_flow_data(flow_data_dir)

    edges = gpd.read_file(edges_file, layer="edges")
    rail_nodes = gpd.read_file(rail_nodes_file, layer="nodes")
    rail_node_ids = rail_nodes[(rail_nodes["asset_type"] == "station") & (rail_nodes["status"] == "Functional")]["node_id"].values.tolist()

    edge_fail_results = []
    for node_fail in rail_node_ids:
        edge_fail: pd.DataFrame = edges[(edges["from_node"] == node_fail) | (edges["to_node"] == node_fail)]
        if len(edge_fail.index) > 0:

            logging.info(f"Failing {node_fail} and adjacent edges")
            # we only remove (extant) edges, but target node_id is included as first entry and becomes label for results row
            to_fail: list[str] = [node_fail] + edge_fail["edge_id"].values.tolist()
            for network in network_dictionary.values():
                edge_fail_results += igraph_scenario_edge_failures_premade_network(
                    network["network"].copy(),
                    to_fail,
                    network["flows"],
                    network["edge_indexes"],
                    "edge_path",
                    "time",
                )

    logging.info("Done failing edges")

    if not edge_fail_results:
        edge_fail_results = pd.DataFrame([], columns=["edge_id", "origin_id", "destination_id", "new_cost", "no_access"])
        logging.info("Node & edge removal had no effect(!)")
    else:
        edge_fail_results = pd.DataFrame(edge_fail_results)
        logging.info(f"Failure results:\n{edge_fail_results}")

    logging.info("Calculating resulting economic losses")
    edge_fail_results.rename(columns={"edge_id": "node_id"}, inplace=True)
    losses = economic_losses_from_network_damage(
        edge_fail_results,
        "node_id",
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
