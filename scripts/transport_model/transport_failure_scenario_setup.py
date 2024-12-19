"""This script allows us to select and parallelise the failure analysis simulations on a server with multipe core processors
    The failure simulations selected here are for:
        Single node and edge failure of the transport networks 
    
    The script first reads in the single node/edge data which we want to fail 
        And then divides the numbers of selected scenarios into paritions which will be selected as initiating failure scenarios
        
        Example output - 0,90
        - First node to sample for initiating failure - The one at location 0 on the sample list
        - Last node to sample for initiating failure - The one at location 90 on the sample list  
        
        We are mainly selecting the first 90 nodes of the road and rail network one-by-one and failing them
        
        All partitions are first written into text files named as parallel_network_scenario_selection.txt 
        Example: parallel_transport_scenario_selection.txt 

        Example output in file: parallel_transport_scenario_selection.txt 
            0,90
            90,181
            181,272
            272,363
            363,453
            453,544
            544,635
        
        Each of these lines is a batch of scenarios that are run on different processors in parallel
"""

import json
import logging
import os
import pandas as pd
import sys

import geopandas as gpd
import numpy as np

from jamaica_infrastructure.transport.utils import load_config
import jamaica_infrastructure.transport.flow as tf


def collate_flow_data(results_path: str, processed_data_path: str):
    """
    Collate flow data and network representations, write to disk in one place.
    This data is used as input to the single link failure analysis.
    """

    logging.info("Reading trade and labour flow data")

    """Get the underlying flow networks
    """
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
            results_path, "flow_mapping", "sector_imports_exports_to_ports_flows.gpkg"
        ),
        layer="edges",
    )
    # print (trade_flow_edges)
    trade_sectors = [
        s
        for s in list(set(trade_flow_edges["from_mode"].values.tolist()))
        if s not in ["road", "rail", "port", "air"]
    ]
    # print (trade_sectors)
    labour_flow_edges = gpd.read_file(
        os.path.join(
            processed_data_path, "networks", "transport", "multi_modal_network.gpkg"
        ),
        layer="edges",
    )
    labour_flow_edges = labour_flow_edges[
        (labour_flow_edges["from_mode"] == "road")
        & (labour_flow_edges["to_mode"] == "road")
    ][columns]
    # print (labour_flow_edges)
    # G = ig.Graph.TupleList(network.itertuples(index=False), edge_attrs=list(network.columns)[2:])

    trade_flows = pd.read_parquet(
        os.path.join(results_path, "flow_mapping", "sector_to_ports_flow_paths.pq")
    )

    labour_flows = pd.read_parquet(
        os.path.join(
            results_path,
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
    labour_flow_path_indexes = tf.get_flow_paths_indexes_of_edges(
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
        sector_flow_path_indexes = tf.get_flow_paths_indexes_of_edges(
            sector_flows, "edge_path"
        )
        network_dictionary["trade"] = {
            "network": sector_network,
            "flows": sector_flows,
            "edge_indexes": pd.Series(sector_flow_path_indexes).to_frame().rename(columns={0: "edge_indexes"}),
        }
        del sector_network, sector_flows, sector_flow_path_indexes

    output_parent_dir = os.path.join(results_path, "transport_failures", "nominal")
    for flow_type, flow_data in network_dictionary.items():
        output_flow_dir = os.path.join(output_parent_dir, flow_type)
        if not os.path.exists(output_flow_dir):
            os.makedirs(output_flow_dir)
        flow_data["network"].to_parquet(os.path.join(output_flow_dir, "network.gpq"))
        flow_data["flows"].to_parquet(os.path.join(output_flow_dir, "flows.pq"))
        flow_data["edge_indexes"].to_parquet(os.path.join(output_flow_dir, "edge_indexes.pq"))

    all_flows.to_parquet(os.path.join(output_parent_dir, "all_flows.pq"))

    with open(os.path.join(output_parent_dir, "trade", "trade_sectors.json"), "w") as fp:
        json.dump(trade_sectors, fp, indent=2)

    return


def main(config, n_chunks):
    processed_data_path = config["paths"]["data"]
    results_path = config["paths"]["output"]

    # collate flow data for labour and trade, write to disk in one location for
    # failure processes to read from
    collate_flow_data(results_path, processed_data_path)

    edges = gpd.read_file(
        os.path.join(
            processed_data_path, "networks", "transport", "multi_modal_network.gpkg"
        ),
        layer="edges",
    )
    rail_edges = edges[(edges["from_mode"] == "rail") & (edges["to_mode"] == "rail")][
        "edge_id"
    ].values.tolist()
    road_edges = edges[(edges["from_mode"] == "road") & (edges["to_mode"] == "road")][
        "edge_id"
    ].values.tolist()
    edge_fail = rail_edges + road_edges

    num_values = np.linspace(0, len(edge_fail) - 1, n_chunks)
    scenarios_path = os.path.join(results_path, "transport_failures", "parallel_transport_scenario_selection.txt")
    with open(scenarios_path, "w+") as f:
        for n in range(len(num_values) - 1):
            f.write(
                "{},{},{}\n".format(
                    "fail edges", int(num_values[n]), int(num_values[n + 1])
                )
            )

    scenarios_resampled_path = os.path.join(
        results_path,
        "transport_failures",
        "parallel_transport_scenario_selection_resample.txt"
    )
    with open(scenarios_resampled_path, "w+") as f:
        with open(scenarios_path) as t:
            for line in t:
                line = line.strip("\n")
                file = os.path.join(
                    results_path,
                    "transport_failures",
                    f"single_link_failures_scenarios_{line.split(',')[1]}_{line.split(',')[2]}.csv",
                )
                print(file)
                if os.path.exists(file) is False:
                    f.write(f"{line}\n")


if __name__ == "__main__":

    _, n_chunks = sys.argv

    logging.basicConfig(format="%(asctime)s %(process)d %(filename)s %(message)s", level=logging.INFO)
    CONFIG = load_config()
    main(CONFIG, int(n_chunks))
