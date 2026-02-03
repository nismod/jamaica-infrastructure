import json
import logging
import os

import click
import geopandas as gpd
import pandas as pd

from jamaica_infrastructure.transport.flow import get_flow_paths_indexes_of_edges


@click.command()
@click.version_option("1.0.0")
@click.option(
    "--results-dir",
    "-r",
    required=True,
    help="Designated results directory.",
)
@click.option(
    "--processed-data-dir",
    "-p",
    required=True,
    help="Designated processed data directory.",
)
def collate_flow_data(results_dir: str, processed_data_dir: str):
    """
    Collate flow data and network representations, write to disk in one place.
    This data is used as input to the single link failure analysis.
    """

    logging.info("Reading trade and labour flow data")

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
    trade_flow_edges = gpd.read_file(
        os.path.join(
            results_dir, "flow_mapping", "sector_imports_exports_to_ports_flows.gpkg"
        ),
        layer="edges",
    )
    trade_sectors = [
        s
        for s in list(set(trade_flow_edges["from_mode"].values.tolist()))
        if s not in ["road", "rail", "port", "air"]
    ]
    labour_flow_edges = gpd.read_file(
        os.path.join(
            processed_data_dir, "networks", "transport", "multi_modal_network.gpkg"
        ),
        layer="edges",
    )
    labour_flow_edges = labour_flow_edges[
        (labour_flow_edges["from_mode"] == "road")
        & (labour_flow_edges["to_mode"] == "road")
    ][columns]

    trade_flows = pd.read_parquet(
        os.path.join(results_dir, "flow_mapping", "sector_to_ports_flow_paths.pq")
    )

    labour_flows = pd.read_parquet(
        os.path.join(
            results_dir,
            "flow_mapping",
            "labour_to_sectors_trips_and_activity.pq",
        )
    )

    logging.info("Combining trade and labour flows")

    all_flows = pd.concat([trade_flows, labour_flows], axis=0, ignore_index=True)
    all_flows = all_flows[
        ["origin_id", "destination_id", "gcost", "working_trips", "GDP_to_trips"]
        + [f"{t}_trade" for t in trade_sectors]
    ]
    labour_flow_path_indexes = get_flow_paths_indexes_of_edges(
        labour_flows, "edge_path"
    )
    network_dictionary = {}
    network_dictionary["labour"] = {
        "network": labour_flow_edges,
        "flows": labour_flows,
        "edge_indexes": pd.Series(labour_flow_path_indexes).to_frame().rename(columns={0: "edge_indexes"}),
    }
    del labour_flow_edges, labour_flows, labour_flow_path_indexes

    for t in trade_sectors:
        logging.info(f"Trade sector: {t}")
        if t != "C":
            sector_network = trade_flow_edges[
                (trade_flow_edges["from_mode"].isin([t, "port", "air", "road"]))
                & (trade_flow_edges["from_mode"].isin([t, "port", "air", "road"]))
            ]
        else:
            sector_network = trade_flow_edges[
                (trade_flow_edges["from_mode"].isin([t, "port", "air", "rail", "road"]))
                & (
                    trade_flow_edges["from_mode"].isin(
                        [t, "port", "air", "rail", "road"]
                    )
                )
            ]
        trade_flows[f"{t}_trade"] = trade_flows[f"{t}_trade"].fillna(0)
        sector_flows = trade_flows[trade_flows[f"{t}_trade"] > 0][
            ["origin_id", "destination_id", "edge_path", "gcost", f"{t}_trade"]
        ].reset_index()
        # TODO: it would be nice if this function just gave us a dataframe we
        # could easily serialise to parquet, rather than a dict
        sector_flow_path_indexes = get_flow_paths_indexes_of_edges(
            sector_flows, "edge_path"
        )
        network_dictionary[f"trade_{t}"] = {
            "network": sector_network,
            "flows": sector_flows,
            "edge_indexes": pd.Series(sector_flow_path_indexes).to_frame().rename(columns={0: "edge_indexes"}),
        }

    logging.info("Writing networks, flows and paths to disk.")

    output_parent_dir = os.path.join(results_dir, "transport_failures", "nominal")
    for flow_type, flow_data in network_dictionary.items():
        output_flow_dir = os.path.join(output_parent_dir, flow_type)
        if not os.path.exists(output_flow_dir):
            os.makedirs(output_flow_dir)
        flow_data["network"].to_parquet(os.path.join(output_flow_dir, "network.gpq"))
        flow_data["flows"].to_parquet(os.path.join(output_flow_dir, "flows.pq"))
        flow_data["edge_indexes"].to_parquet(os.path.join(output_flow_dir, "edge_indexes.pq"))

    all_flows.to_parquet(os.path.join(output_parent_dir, "all_flows.pq"))

    logging.info("Writing trade sectors to disk.")

    with open(os.path.join(output_parent_dir, "trade_sectors.json"), "w") as fp:
        json.dump(trade_sectors, fp, indent=2)

    return


if __name__ == "__main__":

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)
    collate_flow_data()
