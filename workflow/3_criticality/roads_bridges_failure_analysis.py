import logging

import click
import geopandas as gpd
import pandas as pd
from tqdm import tqdm

from jamaica_infrastructure.transport.flow import (
    igraph_scenario_edge_failures_premade_network,
    read_flow_data,
)

tqdm.pandas()
epsg_jamaica = 3448


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
    "--road-nodes-file",
    "-r",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Road nodes GeoPackage file",
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
def main(*, flow_data_dir, edges_file, road_nodes_file, output_path):
    """
    Calculate economic costs of removing road bridge nodes and adjacent edges
    from multi-modal transport network.
    """
    # 0.6 - 2.1% of the value per day
    # so we need to make an assumption on the average wage per working person,
    # say 200 USD per day. Then if a road is disrupted which has 1000 daily trips and
    # they have to be rerouted with an hour, the cost would be: 0.4 * 200 * 1/24 * 100 = 333 USD
    # So corrected for inflation in 2019 values, this would be 1.2-2.9 USD per hour of value of time for business related trips

    hourly_wage = 0.4 * (1 + 0.454) * 235.25  # Between 200 - 500 JMD for 2012 stats, 45.4% inflation in currency
    trade_effect = 0.02  # 2% of the value of trade will be affected by rerouting

    logging.info("Read nominal flow data")
    network_dictionary, all_flows, trade_sectors = read_flow_data(flow_data_dir)

    edges: gpd.GeoDataFrame = gpd.read_file(edges_file, layer="edges")
    road_nodes: gpd.GeoDataFrame = gpd.read_file(road_nodes_file, layer="nodes")
    bridge_node_ids: list[str] = road_nodes[road_nodes["asset_type"] == "bridge"]["node_id"].values.tolist()

    edge_fail_results = []
    for bridge_fail in bridge_node_ids:
        # we chose to make the from_node _only_ of a bridge edge a bridge node
        edge_fail = edges[edges["from_node"] == bridge_fail]
        if len(edge_fail.index) > 0:

            logging.info(f"Failing {bridge_fail} and adjacent edges")
            bridge_edge_ids = edge_fail["edge_id"].values.tolist()
            for networks in network_dictionary.values():
                edge_fail_results += igraph_scenario_edge_failures_premade_network(
                    networks["network"].copy(),
                    bridge_edge_ids,
                    networks["flows"],
                    networks["edge_indexes"],
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
    edge_fail_results = pd.merge(edge_fail_results, all_flows, how="left", on=["origin_id", "destination_id"]).fillna(0)
    edge_fail_results["total_trade"] = edge_fail_results[[f"{t}_trade" for t in trade_sectors]].sum(axis=1)

    edge_fail_results["time_loss"] = (1 - edge_fail_results["no_access"]) * (edge_fail_results["new_cost"] - edge_fail_results["gcost"])
    edge_fail_results["labour_rerouting_loss"] = hourly_wage * edge_fail_results["time_loss"] * edge_fail_results["working_trips"]
    edge_fail_results["trade_rerouting_loss"] = trade_effect * edge_fail_results["time_loss"] * edge_fail_results["total_trade"]
    edge_fail_results["labour_gdp_loss"] = edge_fail_results["no_access"] * edge_fail_results["GDP_to_trips"]
    edge_fail_results["trade_loss"] = edge_fail_results["no_access"] * edge_fail_results["total_trade"]

    index_cols = ["node_id", "no_access"]
    losses = (
        edge_fail_results.loc[
            :,
            index_cols
            + [
                "time_loss",
                "labour_rerouting_loss",
                "trade_rerouting_loss",
                "labour_gdp_loss",
                "trade_loss",
            ],
        ]
        .groupby(index_cols)
        .sum()
        .reset_index()
    )

    rerouting_times_min = edge_fail_results.loc[:, index_cols + ["time_loss"]].groupby(index_cols).min().reset_index()
    rerouting_times_min = rerouting_times_min.rename(columns={"time_loss": "min_trip_time_loss"})
    rerouting_times_max = edge_fail_results.loc[:, index_cols + ["time_loss"]].groupby(index_cols).max().reset_index()
    rerouting_times_max = rerouting_times_max.rename(columns={"time_loss": "max_trip_time_loss"})
    rerouting_times_mean = edge_fail_results.loc[:, index_cols + ["time_loss"]].groupby(index_cols).mean().reset_index()
    rerouting_times_mean = rerouting_times_mean.rename(columns={"time_loss": "mean_trip_time_loss"})

    losses = pd.merge(losses, rerouting_times_min.drop(columns=["no_access"]), how="left", on=["node_id"])
    losses = pd.merge(losses, rerouting_times_max.drop(columns=["no_access"]), how="left", on=["node_id"])
    losses = pd.merge(losses, rerouting_times_mean.drop(columns=["no_access"]), how="left", on=["node_id"])

    logging.info("Writing results to disk")
    losses.to_csv(output_path, index=False)


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)
    main()
